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
    NoReconnections

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

- `ψs`: current streamfunction at vortex filament nodes;

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
    ψs    :: Velocities  # not really velocity, but streamfunction
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
    rhs!     :: RHSFunction     # function for estimating filament velocities (and sometimes streamfunction) from their positions
end

get_dt(iter::VortexFilamentSolver) = iter.time.dt
get_t(iter::VortexFilamentSolver) = iter.time.t

"""
    init(prob::VortexFilamentProblem, scheme::TemporalScheme; dt::Real, kws...) -> VortexFilamentSolver

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
        prob::VortexFilamentProblem, scheme::TemporalScheme;
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
    ψs = similar(vs)
    fs_sol = alias_u0 ? fs : copy.(fs)
    cache_bs = BiotSavart.init_cache(prob.p, fs_sol; timer)
    cache_timestepper = init_cache(scheme, fs, vs)

    # Wrap functions with the timer, so that timings are estimated each time the function is called.
    advect! = timer(advect_filaments!, "advect_filaments!")
    rhs! = timer(update_values_at_nodes!, "update_values_at_nodes!")
    callback_ = timer(callback, "callback")

    if adaptivity !== NoAdaptivity() && !can_change_dt(scheme)
        throw(ArgumentError(lazy"temporal scheme $scheme doesn't support adaptibility; set `adaptivity = NoAdaptivity()` or choose a different scheme"))
    end

    time = TimeInfo(nstep = 0, t = first(tspan), dt = dt, dt_prev = dt)

    iter = VortexFilamentSolver(
        prob, fs_sol, vs, ψs, time, dtmin, refinement, adaptivity, reconnect,
        cache_bs, cache_timestepper, callback_, timer,
        advect!, rhs!,
    )

    finalise_step!(iter)

    iter
end

function _update_values_at_nodes!(
        ::Val{:full},
        fields::NamedTuple{Names, NTuple{N, V}},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    ) where {Names, N, V <: VectorOfVectors}
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(true))
    fields
end

# Compute slow component only.
# This is generally called in IMEX-RK substeps, where only the velocity (and not the
# streamfunction) is needed.
# We assume that the "slow" component is everything but LIA term when evolving the
# Biot-Savart law.
# This component will be treated explicitly by IMEX schemes.
function _update_values_at_nodes!(
        ::Val{:slow},
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(false))
    fields
end

# Compute fast component only.
# This is generally called in IMEX-RK substeps, where only the velocity (and not the
# streamfunction) is needed.
# We assume that the "fast" component is the LIA term when evolving the Biot-Savart law.
# This component will be treated implicitly by IMEX schemes.
function _update_values_at_nodes!(
        ::Val{:fast},
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    vs_all = fields.velocity
    (; prob, to,) = iter
    params = prob.p  # Biot-Savart parameters
    (; Γ, a, Δ,) = params.common
    # Note: we must use the same quadrature as used when using component = Val(:full)
    (; quad,) = params.shortrange
    prefactor = Γ / (4π)
    @timeit to "LIA term (only)" begin
        @inbounds for (f, vs) ∈ zip(fs, vs_all)
            for i ∈ eachindex(f, vs)
                vs[i] = BiotSavart.local_self_induced(
                    BiotSavart.Velocity(), f, i, prefactor;
                    a, Δ, quad,
                )
            end
        end
    end
    fields
end

# This is the most general variant which should be called by timesteppers.
# For now we don't use the time, but we might in the future...
function update_values_at_nodes!(
        fields::NamedTuple, fs, t::Real, iter;
        component = Val(:full),  # compute slow + fast components by default
    )
    _update_values_at_nodes!(component, fields, fs, iter)
end

# Case where only the velocity is passed (generally used in internal RK substeps).
function update_values_at_nodes!(vs::VectorOfVectors, args...; kws...)
    fields = (velocity = vs,)
    update_values_at_nodes!(fields, args...; kws...)
    vs
end

# Here buf is usually a VectorOfFilaments or a vector of vectors of velocities
# (containing the velocities of all filaments).
# Resizes the higher-level vector without resizing the individual vectors it contains.
function resize_container!(buf, fs::VectorOfFilaments)
    i = lastindex(buf)
    N = lastindex(fs)
    i === N && return buf
    while i < N
        i += 1
        push!(buf, similar(first(buf), length(fs[i])))
    end
    while i > N
        i -= 1
        pop!(buf)
    end
    @assert length(fs) == length(buf)
    buf
end

function resize_contained_vectors!(buf, fs::VectorOfFilaments)
    for (i, f) ∈ pairs(fs)
        N = length(f)
        resize!(buf[i], N)
    end
    buf
end

function refine!(f::AbstractFilament, refinement::RefinementCriterion)
    nref = Filaments.refine!(f, refinement)  # we assume update_coefficients! was already called
    n = 0
    nmax = 1  # maximum number of extra refinement iterations (to avoid infinite loop...)
    while nref !== (0, 0) && n < nmax
        n += 1
        @debug "Added/removed $nref nodes"
        if !Filaments.check_nodes(Bool, f)
            return nothing  # filament should be removed (too small / not enough nodes)
        end
        nref = Filaments.refine!(f, refinement)
    end
    n
end

function _advect_filament!(f::T, fbase::T, vs, dt) where {T <: AbstractFilament}
    Xs = nodes(f)
    Ys = nodes(fbase)
    for i ∈ eachindex(Xs, Ys, vs)
        @inbounds Xs[i] = Ys[i] + dt * vs[i]
    end
    # Compute interpolation coefficients making sure that knots are preserved (and not
    # recomputed from new locations).
    Filaments.update_coefficients!(f; knots = knots(fbase))
    nothing
end

function _advect_filament!(
        f::AbstractFilament, fbase::Nothing, vs::VectorOfVelocities, dt::Real,
    )
    Xs = nodes(f)
    for i ∈ eachindex(Xs, vs)
        @inbounds Xs[i] = Xs[i] + dt * vs[i]
    end
    Filaments.update_coefficients!(f)  # note: in this case we recompute the knots
    nothing
end

# This should be called right after positions have been updated by the timestepper, but
# before reconnections happen.
# Here `fields` is a tuple containing the fields associated to each filament node.
# Usually this is `(vs, ψs)`, which contain the velocities and streamfunction values at all
# filament nodes.
function after_advection!(fs, fields::Tuple; kws...)
    @inbounds for i ∈ reverse(eachindex(fs))
        fields_i = map(vs -> @inbounds(vs[i]), fields)  # typically this is (vs[i], ψs[i])
        ret = _after_advection!(fs[i], fields_i; kws...)
        if ret === nothing
            # This should only happen if the filament was "refined" (well, unrefined in this case).
            popat!(fs, i)
            map(vs -> popat!(vs, i), fields)
        end
    end
    fs
end

function _after_advection!(f::AbstractFilament, fields::Tuple; L_fold, refinement,)
    if Filaments.fold_periodic!(f, L_fold)
        Filaments.update_coefficients!(f)  # only called if nodes were modified by fold_periodic!
    end
    refinement_steps = refine!(f, refinement)
    if refinement_steps === nothing  # filament should be removed (too small / not enough nodes)
        return nothing
    elseif refinement_steps > 0  # refinement was performed
        map(vs -> resize!(vs, length(f)), fields)
    end
    f
end

# The possible values of fbase are:
# - `nothing` (default), used when the filaments `fs` should be actually advected with the
#   chosen velocity at the end of a timestep;
# - some list of filaments similar to `fs`, used to advect from `fbase` to `fs` with the
#   chosen velocity. This is generally done from within RK stages, and in this case things
#   like filament refinement are not performed.
advect_filaments!(fs, vs, dt; fbase = nothing, kws...) =
    _advect_filaments!(fs, fbase, vs, dt; kws...)

# Variant called at the end of a timestep.
function _advect_filaments!(fs, fbase::Nothing, vs, dt)
    for i ∈ eachindex(fs, vs)
        @inbounds _advect_filament!(fs[i], nothing, vs[i], dt)
    end
    fs
end

# Variant called from within RK stages.
function _advect_filaments!(fs::T, fbase::T, vs, dt) where {T}
    for i ∈ eachindex(fs, fbase, vs)
        @inbounds _advect_filament!(fs[i], fbase[i], vs[i], dt)
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
    (; fs, vs, ψs, prob, time, dtmin, refinement, advect!, rhs!, to,) = iter
    t_end = iter.prob.tspan[2]
    if time.dt < dtmin && time.t + time.dt < t_end
        error(lazy"current timestep is too small ($(time.dt) < $(dtmin)). Stopping.")
    end
    # Note: the timesteppers assume that iter.vs already contains the velocity
    # induced by the filaments at the current timestep.
    update_velocities!(
        rhs!, advect!, iter.cache_timestepper, iter,
    )
    L_fold = periods(prob.p)  # box size (periodicity)
    advect!(fs, vs, time.dt)
    after_advection!(fs, (vs, ψs); L_fold, refinement)
    time.t += time.dt
    time.nstep += 1
    finalise_step!(iter)
    iter
end

# Here fields is usually (vs, ψs)
@inline function reconnect_callback((fs, fields), f, i, mode::Symbol)
    if mode === :removed
        @debug lazy"Filament was removed at index $i"
        # @assert i ≤ lastindex(fields[1]) == lastindex(fs) + 1
        map(vs -> popat!(vs, i), fields)
    elseif mode === :appended
        @debug lazy"Filament was appended at index $i"
        @assert f === fs[i]
        # @assert i == lastindex(fields[1]) + 1
        map(vs -> push!(vs, similar(first(vs), length(f))), fields)
    elseif mode === :modified
        @debug lazy"Filament was modified at index $i"
        @assert f === fs[i]
        # @assert i ≤ lastindex(fs) == lastindex(fields[1])
        map(vs -> resize!(vs[i], length(f)), fields)
    end
    nothing
end

function reconnect!(iter::VortexFilamentSolver)
    (; vs, ψs, fs, reconnect,) = iter
    fields = (vs, ψs)
    Filaments.reconnect!(reconnect, fs) do f, i, mode
        reconnect_callback((fs, fields), f, i, mode)
    end :: Int  # returns the number of reconnections
end

# Called whenever filament positions have just been initialised or updated.
function finalise_step!(iter::VortexFilamentSolver)
    (; vs, ψs, fs, time, callback, adaptivity, rhs!, to,) = iter

    # Perform reconnections, possibly changing the number of filaments.
    @timeit to "reconnect!" number_of_reconnections = reconnect!(iter)
    @debug lazy"Number of reconnections: $number_of_reconnections"

    @assert eachindex(fs) == eachindex(vs)

    # Check if there are filaments to be removed (typically with ≤ 3 discretisation points,
    # but this depends on the precise discretisation method). Note that filaments may also
    # be removed during reconnections for the same reason.
    for i ∈ reverse(eachindex(fs))
        if !Filaments.check_nodes(Bool, fs[i])
            # Remove filament and its associated vector of velocities
            popat!(fs, i)
            popat!(vs, i)
        end
    end

    isempty(fs) && error("all vortices disappeared!")  # TODO nicer way to handle this?

    # Update velocities and streamfunctions to the next timestep (and first RK step).
    # Note that we only compute the streamfunction at full steps, and not in the middle of
    # RK substeps.
    # TODO make computation of ψ optional?
    fields = (velocity = vs, streamfunction = ψs,)

    # Note: here we always include the LIA terms, even when using IMEX schemes.
    # This must be taken into account by IMEX scheme implementations.
    rhs!(fields, fs, time.t, iter; component = Val(:full))

    # This is mainly useful for visualisation (and it's quite cheap).
    for u ∈ vs
        Filaments.pad_periodic!(u)
    end
    for u ∈ ψs
        Filaments.pad_periodic!(u)
    end

    callback(iter)
    time.dt_prev = time.dt
    time.dt = estimate_timestep(adaptivity, iter)  # estimate dt for next timestep

    iter
end

end
