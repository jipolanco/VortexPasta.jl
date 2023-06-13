"""
    Timestepping

Module defining timestepping solvers for vortex filament simulations.
"""
module Timestepping

export init, solve!, step!, VortexFilamentProblem

using ..BasicTypes: Vec3, VectorOfVectors

using ..Filaments:
    Filaments,
    AbstractFilament,
    nodes,
    segments,
    knots,

    RefinementCriterion,
    NoRefinement,

    ReconnectionCriterion,
    NoReconnections,

    update_coefficients_before_refinement,
    update_coefficients_after_refinement

using ..BiotSavart:
    BiotSavart,
    ParamsBiotSavart,
    BiotSavartCache,
    VectorOfFilaments,
    VectorOfPositions,
    VectorOfVelocities,
    AllFilamentVelocities,
    periods

# Reuse same init, solve! and step! functions from the SciML ecosystem, to avoid clashes.
# See https://docs.sciml.ai/CommonSolve/stable/
import CommonSolve: init, solve!, step!

using TimerOutputs: TimerOutput, @timeit

abstract type AbstractProblem end
abstract type AbstractSolver end

include("timesteppers/timesteppers.jl")
include("adaptivity.jl")

"""
    VortexFilamentProblem
    VortexFilamentProblem(fs::AbstractVector{<:AbstractFilament}, tspan, p::ParamsBiotSavart)

Define a vortex filament problem.

Arguments:

- `fs`: initial vortex positions;

- `tspan = (t_begin, t_end)`: time span;

- `p`: Biot–Savart parameters (see [`ParamsBiotSavart`](@ref)).

See [`init`](@ref) for initialising a solver from a `VortexFilamentProblem`.
"""
struct VortexFilamentProblem{
        Filaments <: VectorOfFilaments,
        Params <: ParamsBiotSavart,
    } <: AbstractProblem
    fs    :: Filaments
    tspan :: NTuple{2, Float64}
    p     :: Params
end

"""
    TimeInfo

Contains information on the current time and timestep of a solver.

Some useful fields are:

- `t`: current time;

- `dt`: timestep to be used in next iteration;

- `dt_prev` : timestep used in the last performed iteration;

- `nstep`: number of timesteps performed until now.
"""
@kwdef mutable struct TimeInfo
    nstep   :: Int
    t       :: Float64
    dt      :: Float64
    dt_prev :: Float64
end

"""
    VortexFilamentSolver

Contains the instantaneous state of a vortex filament simulation.

Must be constructed using [`init`](@ref).

Some useful fields are:

- `prob`: associated [`VortexFilamentProblem`](@ref) (including Biot–Savart parameters);

- `fs`: current state of vortex filaments in the system;

- `vs`: current velocity of vortex filament nodes;

- `time`: a [`TimeInfo`](@ref) object containing information such as the current time and timestep;

- `to`: a `TimerOutput`, which records the time spent on different functions;

- `cache_bs`: the Biot–Savart cache, which contains data from short- and
  long-range computations.
"""
struct VortexFilamentSolver{
        Problem <: VortexFilamentProblem,
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfVectors{<:Vec3},
        Refinement <: RefinementCriterion,
        Adaptivity <: AdaptivityCriterion,
        Reconnect <: ReconnectionCriterion,
        CacheBS <: BiotSavartCache,
        CacheTimestepper <: TemporalSchemeCache,
        Timer <: TimerOutput,
        Callback <: Function,
        AdvectFunction <: Function,
        RHSFunction <: Function,
    } <: AbstractSolver
    prob  :: Problem
    fs    :: Filaments
    vs    :: Velocities
    time  :: TimeInfo
    dtmin :: Float64
    refinement        :: Refinement
    adaptivity        :: Adaptivity
    reconnect         :: Reconnect
    cache_bs          :: CacheBS
    cache_timestepper :: CacheTimestepper
    callback :: Callback
    to       :: Timer
    advect!  :: AdvectFunction  # function for advecting filaments with a known velocity
    rhs!     :: RHSFunction     # function for estimating filament velocities from their positions
end

get_dt(iter::VortexFilamentSolver) = iter.time.dt
get_t(iter::VortexFilamentSolver) = iter.time.t

"""
    init(prob::VortexFilamentProblem, scheme::ExplicitTemporalScheme; dt::Real, kws...) -> VortexFilamentSolver

Initialise vortex filament problem.

Returns a [`VortexFilamentSolver`](@ref) which can be advanced in time using
either [`step!`](@ref) or [`solve!`](@ref).

# Mandatory keyword arguments

- `dt::Real`: the simulation timestep.

# Optional keyword arguments

- `alias_u0 = true`: if `true` (default), the solver is allowed to modify the
  initial vortex filaments.

- `refinement = NoRefinement()`: criterion used for adaptive refinement of vortex
  filaments. See [`RefinementCriterion`](@ref) for a list of possible criteria.

- `reconnect = NoReconnections()`: criterion used to perform vortex reconnections.
  See [`ReconnectionCriterion`](@ref) for a list of possible criteria.

- `adaptivity = NoAdaptivity()`: criterion used for adaptively setting the timestep `dt`.
  See [`AdaptivityCriterion`](@ref) for a list of possible criteria.

- `dtmin = 0.0`: minimum `dt` for adaptive timestepping. If `dt < dtmin`, the
  solver is stopped with an error.

- `callback`: a function to be called at the end of each timestep. The function
  must accept a single argument `iter::VortexFilamentSolver`.

- `timer = TimerOutput("VortexFilament")`: an optional `TimerOutput` for
  recording the time spent on different functions.
"""
function init(
        prob::VortexFilamentProblem, scheme::ExplicitTemporalScheme;
        alias_u0 = true,   # same as in OrdinaryDiffEq.jl
        dt::Real,
        dtmin::Real = 0.0,
        refinement::RefinementCriterion = NoRefinement(),
        reconnect::ReconnectionCriterion = NoReconnections(),
        adaptivity::AdaptivityCriterion = NoAdaptivity(),
        callback::F = identity,
        timer = TimerOutput("VortexFilament"),
    ) where {F <: Function}
    (; fs, tspan,) = prob
    vs_data = [similar(nodes(f)) for f ∈ fs] :: AllFilamentVelocities
    vs = VectorOfVectors(vs_data)
    fs_sol = alias_u0 ? fs : copy.(fs)
    cache_bs = BiotSavart.init_cache(prob.p; timer)
    cache_timestepper = init_cache(scheme, fs, vs)

    # Wrap functions with the timer, so that timings are estimated each time the function is called.
    advect! = timer(advect_filaments!, "advect_filaments!")
    rhs! = timer(vortex_velocities!, "vortex_velocities!")
    callback_ = timer(callback, "callback")

    if adaptivity !== NoAdaptivity() && !can_change_dt(scheme)
        throw(ArgumentError(lazy"temporal scheme $scheme doesn't support adaptibility; set `adaptivity = NoAdaptivity()` or choose a different scheme"))
    end

    time = TimeInfo(nstep = 0, t = first(tspan), dt = dt, dt_prev = dt)

    iter = VortexFilamentSolver(
        prob, fs_sol, vs, time, dtmin, refinement, adaptivity, reconnect,
        cache_bs, cache_timestepper, callback_, timer,
        advect!, rhs!,
    )

    after_advection!(iter)

    iter
end

function _vortex_velocities!(
        vs::VectorOfVectors,
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    (; cache_bs,) = iter
    BiotSavart.velocity_on_nodes!(vs.u, cache_bs, fs)
    vs
end

# This is the most general variant which should be called by timesteppers.
# For now we don't use the time, but we might in the future...
vortex_velocities!(vs, fs, t::Real, iter) = _vortex_velocities!(vs, fs, iter)

function refine!(f::AbstractFilament, refinement::RefinementCriterion)
    nref = Filaments.refine!(f, refinement)  # we assume update_coefficients! was already called
    n = 0
    while nref !== (0, 0)
        n += 1
        @debug "Added/removed $nref nodes"
        if !Filaments.check_nodes(Bool, f)
            return nothing  # filament should be removed (too small / not enough nodes)
        end
        if update_coefficients_after_refinement(f)
            Filaments.update_coefficients!(f)
        end
        nref = Filaments.refine!(f, refinement)
    end
    n
end

function _advect_filament!(
        f::AbstractFilament, vs::VectorOfVelocities, dt::Real;
        fbase = f, L_fold = nothing, refinement = nothing,
    )
    Xs = nodes(f)
    @assert Filaments.check_nodes(Bool, f)  # filament has enough nodes
    Xs_base = nodes(fbase)
    @assert eachindex(Xs) == eachindex(Xs_base) == eachindex(vs)
    @inbounds for i ∈ eachindex(Xs, Xs_base, vs)
        Xs[i] = Xs_base[i] + dt * vs[i]
    end
    if L_fold !== nothing
        Filaments.fold_periodic!(f, L_fold)
    end
    need_to_update_coefs = true
    if refinement === nothing
        Filaments.update_coefficients!(f)
    else
        # Check whether the filament type requires coefficients to be updated before refining.
        # The answer may be different depending on whether we're using spline or finite difference
        # discretisations for filaments.
        if update_coefficients_before_refinement(f)
            Filaments.update_coefficients!(f)
        end
        refinement_steps = refine!(f, refinement)
        if refinement_steps === nothing  # filament should be removed (too small / not enough nodes)
            return nothing
        elseif refinement_steps > 0  # refinement was performed
            resize!(vs, length(f))
        end
        need_to_update_coefs = false  # already done in `refine!`
    end
    f
end

function advect_filaments!(fs, vs, dt; fbase = fs, kws...)
    for i ∈ reverse(eachindex(fs, fbase, vs))
        ret = _advect_filament!(fs[i], vs[i], dt; fbase = fbase[i], kws...)
        if ret === nothing
            # This should only happen if the filament was "refined" (well, unrefined in this case).
            # This only happens after a full timestep (i.e. not in intermediate RK stages),
            # in which case `fbase` and `fs` are the same.
            @assert fs === fbase
            popat!(fs, i)
            popat!(vs, i)
        end
    end
    fs
end

"""
    solve!(iter::VortexFilamentSolver)

Advance vortex filament solver to the ending time.

See also [`step!`](@ref) for advancing one step at a time.
"""
function solve!(iter::VortexFilamentSolver)
    (; time,) = iter
    t_end = iter.prob.tspan[2]
    while time.t < t_end
        if can_change_dt(iter.cache_timestepper)
            # Try to finish exactly at t = t_end.
            time.dt = min(time.dt, t_end - time.t)
        end
        step!(iter)
    end
    iter
end

"""
    step!(iter::VortexFilamentSolver)

Advance solver by a single timestep.
"""
function step!(iter::VortexFilamentSolver)
    (; fs, vs, prob, time, dtmin, refinement, advect!, rhs!, to,) = iter
    time.dt ≥ dtmin || error(lazy"current timestep is too small ($(time.dt) < $(dtmin)). Stopping.")
    # Note: the timesteppers assume that iter.vs already contains the velocity
    # induced by the filaments at the current timestep.
    update_velocities!(
        rhs!, advect!, iter.cache_timestepper, iter,
    )
    L_fold = periods(prob.p)  # box size (periodicity)
    advect!(fs, vs, time.dt; L_fold, refinement)
    time.t += time.dt
    time.nstep += 1
    after_advection!(iter)
    iter
end

function reconnect!(iter::VortexFilamentSolver)
    (; vs, fs, reconnect,) = iter
    Filaments.reconnect!(reconnect, fs) do f, i, mode
        if mode === :removed
            @debug lazy"Filament was removed at index $i"
            @assert i ≤ lastindex(vs) == lastindex(fs) + 1
            popat!(vs, i)
        elseif mode === :appended
            @debug lazy"Filament was appended at index $i"
            @assert f === fs[i]
            @assert i == lastindex(vs) + 1
            push!(vs, similar(first(vs), length(f)))
        elseif mode === :modified
            @debug lazy"Filament was modified at index $i"
            @assert f === fs[i]
            @assert i ≤ lastindex(fs) == lastindex(vs)
            resize!(vs[i], length(f))
        end
    end
    iter
end

# Called whenever filament positions have just been initialised or updated.
function after_advection!(iter::VortexFilamentSolver)
    (; vs, fs, time, callback, adaptivity, rhs!, to,) = iter

    # Perform reconnections, possibly changing the number of filaments.
    @timeit to "reconnect!" reconnect!(iter)
    @assert length(fs) == length(vs)

    rhs!(vs, fs, time.t, iter)  # update velocities to the next timestep (and first RK step)
    callback(iter)
    time.dt_prev = time.dt
    time.dt = estimate_timestep(adaptivity, iter)  # estimate dt for next timestep

    iter
end

end
