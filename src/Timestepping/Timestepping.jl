"""
    Timestepping

Module defining timestepping solvers for vortex filament simulations.
"""
module Timestepping

export init, solve!, step!, VortexFilamentProblem

using ..BasicTypes: Vec3

using ..Filaments:
    Filaments,
    AbstractFilament,
    nodes,
    segments,
    knots,
    RefinementCriterion,
    BasedOnCurvature

using ..BiotSavart:
    BiotSavart,
    ParamsBiotSavart,
    BiotSavartCache,
    VectorOfFilaments,
    VectorOfPositions,
    VectorOfVelocities,
    AllFilamentVelocities,
    periods

# Reuse same init and solve! functions from the SciML ecosystem, to avoid clashes.
# See https://docs.sciml.ai/CommonSolve/stable/
import CommonSolve: init, solve!
using RecursiveArrayTools: VectorOfArray  # for convenience, to deal with filament velocities

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

- `dt`: timestep for next iteration;

- `to`: a `TimerOutput`, which records the time spent on different functions;

- `cache_bs`: the Biot–Savart cache, which contains data from short- and
  long-range computations.
"""
struct VortexFilamentSolver{
        Problem <: VortexFilamentProblem,
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfArray{<:Vec3},
        Refinement <: RefinementCriterion,
        Adaptivity <: AdaptivityCriterion,
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

- `refinement = BasedOnCurvature(0.35; ℓ_max = 1.0)`: method used for adaptive
  refinement of vortex filaments. See [`BasedOnCurvature`](@ref) for details.

- `adaptivity = NoAdaptivity()`: method used for adaptively setting the timestep `dt`.
  See [`AdaptivityCriterion`](@ref) for a list of possible methods and
  [`BasedOnSegmentLength`](@ref) for one of these methods.

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
        refinement::RefinementCriterion = BasedOnCurvature(0.35; ℓ_max = 1.0),
        adaptivity::AdaptivityCriterion = NoAdaptivity(),
        callback::F = identity,
        timer = TimerOutput("VortexFilament"),
    ) where {F <: Function}
    (; fs, tspan,) = prob
    vs_data = [similar(nodes(f)) for f ∈ fs] :: AllFilamentVelocities
    vs = VectorOfArray(vs_data)
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
        prob, fs_sol, vs, time, dtmin, refinement, adaptivity,
        cache_bs, cache_timestepper, callback_, timer,
        advect!, rhs!,
    )

    after_advection!(iter)

    iter
end

function _vortex_velocities!(
        vs::VectorOfArray,
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

function _advect_filament!(
        f::AbstractFilament, vs::VectorOfVelocities, dt::Real;
        fbase = f, L_fold = nothing, refinement = nothing,
    )
    Xs = nodes(f)
    Xs_base = nodes(fbase)
    @assert eachindex(Xs) == eachindex(Xs_base) == eachindex(vs)
    @inbounds for i ∈ eachindex(Xs)
        Xs[i] = Xs_base[i] + dt * vs[i]
    end
    if L_fold !== nothing
        Filaments.fold_periodic!(f, L_fold)
    end
    if refinement !== nothing
        nref = Filaments.refine!(f, refinement)
        if nref !== (0, 0)
            # @info "Added/removed $nref nodes"
            if sum(nref) ≠ 0
                resize!(vs, length(f))
            end
        end
    end
    Filaments.update_coefficients!(f)
    f
end

function advect_filaments!(fs, vs, dt; fbase = fs, kws...)
    for (f, f₀, v) ∈ zip(fs, fbase, vs)
        _advect_filament!(f, v, dt; fbase = f₀, kws...)
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

# Called whenever filament positions have just been initialised or updated.
function after_advection!(iter::VortexFilamentSolver)
    (; vs, fs, time, callback, adaptivity, rhs!,) = iter
    rhs!(vs, fs, time.t, iter)  # update velocities to the next timestep (and first RK step)
    callback(iter)
    time.dt_prev = time.dt
    time.dt = estimate_timestep(adaptivity, iter)  # estimate dt for next timestep
    iter
end

end
