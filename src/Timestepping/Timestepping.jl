"""
    Timestepping

Module defining timestepping solvers for vortex filament simulations.
"""
module Timestepping

export init, solve!, step!, VortexFilamentProblem,
       ShortRangeTerm, LocalTerm,
       ParamsBiotSavart,                           # from ..BiotSavart
       NoReconnections, ReconnectBasedOnDistance,  # from ..Reconnections
       reset_timer!  # from TimerOutputs

using ..BasicTypes: Vec3, VectorOfVectors

using ..Filaments:
    Filaments,
    AbstractFilament,
    nodes,
    segments,
    knots,
    RefinementCriterion,
    NoRefinement

using ..Reconnections:
    Reconnections,
    AbstractReconnectionCache,
    ReconnectionCriterion,
    NoReconnections,
    ReconnectBasedOnDistance,
    reconnect!

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

using ForwardDiff: ForwardDiff
using TimerOutputs: TimerOutputs, TimerOutput, @timeit, reset_timer!

abstract type AbstractProblem end
abstract type AbstractSolver end

include("timesteppers/timesteppers.jl")
include("adaptivity.jl")
include("problem.jl")

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

function Base.show(io::IO, time::TimeInfo)
    (; nstep, t, dt, dt_prev,) = time
    print(io, "TimeInfo:")
    print(io, "\n - nstep   = ", nstep)
    print(io, "\n - t       = ", t)
    print(io, "\n - dt      = ", dt)
    print(io, "\n - dt_prev = ", dt_prev)
end

function Base.summary(io::IO, time::TimeInfo)
    (; nstep, t, dt, dt_prev,) = time
    print(io, "TimeInfo(nstep = $nstep, t = $t, dt = $dt, dt_prev = $dt_prev)")
end

abstract type FastBiotSavartTerm end

"""
    ShortRangeTerm <: FastBiotSavartTerm

Identifies fast dynamics with the **short-range component** of Ewald splitting.

This is only useful for split timestepping schemes like IMEX or multirate methods.
"""
struct ShortRangeTerm <: FastBiotSavartTerm end

"""
    LocalTerm <: FastBiotSavartTerm

Identifies fast dynamics with the **local (LIA) term** associated to the desingularisation
of the Biot–Savart integral.

This is useful for split timestepping schemes like IMEX or multirate methods.
"""
struct LocalTerm <: FastBiotSavartTerm end

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
        T,
        Problem <: VortexFilamentProblem{T},
        Filaments <: VectorOfVectors{Vec3{T}, <:AbstractFilament{T}},
        Velocities <: VectorOfVectors{Vec3{T}},
        Refinement <: RefinementCriterion,
        Adaptivity <: AdaptivityCriterion,
        CacheReconnect <: AbstractReconnectionCache,
        CacheBS <: BiotSavartCache,
        CacheTimestepper <: TemporalSchemeCache,
        FastTerm <: FastBiotSavartTerm,
        Timer <: TimerOutput,
        Callback <: Function,
        ExternalForcing <: NamedTuple,
        AdvectFunction <: Function,
        RHSFunction <: Function,
    } <: AbstractSolver
    prob  :: Problem
    fs    :: Filaments
    vs    :: Velocities
    ψs    :: Velocities  # not really velocity, but streamfunction
    time  :: TimeInfo
    dtmin :: T
    refinement        :: Refinement
    adaptivity        :: Adaptivity
    reconnect         :: CacheReconnect
    cache_bs          :: CacheBS
    cache_timestepper :: CacheTimestepper
    fast_term         :: FastTerm
    callback :: Callback
    external_forcing :: ExternalForcing  # velocity and streamfunction forcing
    to       :: Timer
    advect!  :: AdvectFunction  # function for advecting filaments with a known velocity
    rhs!     :: RHSFunction     # function for estimating filament velocities (and sometimes streamfunction) from their positions
end

function Base.show(io::IO, iter::VortexFilamentSolver)
    print(io, "VortexFilamentSolver with fields:")
    print(io, "\n - `prob`: ", nameof(typeof(iter.prob)))
    _print_summary(io, iter.fs; pre = "\n - `fs`: ", post =  " -- vortex filaments")
    _print_summary(io, iter.vs; pre = "\n - `vs`: ", post = " -- velocity at nodes")
    _print_summary(io, iter.ψs; pre = "\n - `ψs`: ", post = " -- streamfunction at nodes")
    print(io, "\n - `time`: ")
    summary(io, iter.time)
    print(io, "\n - `dtmin`: ", iter.dtmin)
    print(io, "\n - `refinement`: ", iter.refinement)
    print(io, "\n - `adaptivity`: ", iter.adaptivity)
    print(io, "\n - `reconnect`: ", Reconnections.criterion(iter.reconnect))
    print(io, "\n - `cache_bs`: ")
    summary(io, iter.cache_bs)
    print(io, "\n - `cache_timestepper`: ")
    summary(io, iter.cache_timestepper)
    print(io, "\n - `fast_term`: ", iter.fast_term)
    print(io, "\n - `callback`: Function (`", _printable_function(iter.callback), "`)")
    for (name, func) ∈ pairs(iter.external_forcing)
        func === nothing || print(io, "\n - `external_forcing.$name`: Function (`$func`)")
    end
    print(io, "\n - `advect!`: Function")
    print(io, "\n - `rhs!`: Function")
    print(io, "\n - `to`: ")
    summary(io, iter.to)
end

# This is for convenience: allows doing e.g. `iter.t` instead of `iter.time.t` to get the
# current time.
@inline function Base.getproperty(iter::VortexFilamentSolver, name::Symbol)
    time = getfield(iter, :time)
    if hasproperty(time, name)
        getproperty(time, name)
    else
        getfield(iter, name)
    end
end

function Base.propertynames(iter::VortexFilamentSolver, private::Bool = false)
    (fieldnames(typeof(iter))..., propertynames(iter.time, private)...)
end

_printable_function(f::TimerOutputs.InstrumentedFunction) = f.func

get_dt(iter::VortexFilamentSolver) = iter.time.dt
get_t(iter::VortexFilamentSolver) = iter.time.t

@doc raw"""
    init(prob::VortexFilamentProblem, scheme::TemporalScheme; dt::Real, kws...) -> VortexFilamentSolver

Initialise vortex filament problem.

Returns a [`VortexFilamentSolver`](@ref) which can be advanced in time using
either [`step!`](@ref) or [`solve!`](@ref).

# Mandatory keyword arguments

- `dt::Real`: the simulation timestep.

# Optional keyword arguments

- `alias_u0 = false`: if `true`, the solver is allowed to modify the initial vortex
  filaments.

- `refinement = NoRefinement()`: criterion used for adaptive refinement of vortex
  filaments. See [`RefinementCriterion`](@ref) for a list of possible criteria.

- `reconnect = NoReconnections()`: criterion used to perform vortex reconnections.
  See [`ReconnectionCriterion`](@ref) for a list of possible criteria.

- `adaptivity = NoAdaptivity()`: criterion used for adaptively setting the timestep `dt`.
  See [`AdaptivityCriterion`](@ref) for a list of possible criteria.

- `dtmin = 0.0`: minimum `dt` for adaptive timestepping. If `dt < dtmin`, the
  solver is stopped with an error.

- `fast_term = LocalTerm()`: for IMEX and multirate schemes, this determines what is meant by
  "fast" and "slow" dynamics. This can either be [`LocalTerm`](@ref) or [`ShortRangeTerm`](@ref).
  Note that the default may change in the future!

- `callback`: a function to be called at the end of each timestep. The function
  must accept a single argument `iter::VortexFilamentSolver`.

- `timer = TimerOutput("VortexFilament")`: an optional `TimerOutput` for
  recording the time spent on different functions.

# Extended help

## Adding an external velocity

One can set the **`external_velocity`** keyword argument to impose an external velocity field
``\bm{v}_{\text{f}}(\bm{x}, t)``.
In this case, the total velocity at a point ``\bm{x}`` (which can be on a vortex filament) will be

```math
\bm{v}(\bm{x}) = \bm{v}_{\text{BS}}(\bm{x}) + \bm{v}_{\text{f}}(\bm{x})
```

where ``\bm{v}_{\text{BS}}`` is the velocity obtained from the Biot–Savart law.

This can be a way of applying an external forcing (and injecting energy) to the vortex system.

The external velocity should be given as a function, which should have the signature `v_ext(x⃗::Vec3, t::Real) -> Vec3`.
The function should return a `Vec3` with the velocity at the location `x⃗` and time `t`.
For example, to impose a constant sinusoidal forcing in the ``x`` direction which varies along ``z``,
the `external_velocity` keyword argument could look like:

    external_velocity = (x⃗, t) -> Vec3(0.1 * sin(2 * x⃗.z), 0, 0)

One usually wants the external velocity to be divergence-free, to preserve the
incompressibility of the flow.

### Notes on energy estimation

When setting an external velocity ``\bm{v}_{\text{f}}(\bm{x}, t)``, one may also want to
impose an external streamfunction field ``\bm{ψ}_{\text{f}}(\bm{x}, t)`` via the **`external_streamfunction`**
keyword argument.
Such a field plays no role in the dynamics, but it allows to include the kinetic energy
associated to the external velocity field when calling [`Diagnostics.kinetic_energy_from_streamfunction`](@ref).
If provided, the external streamfunction should satisfy ``\bm{v}_{\text{f}} = \bm{∇} × \bm{ψ}_{\text{f}}``.

More precisely, when applying an external velocity field, the total kinetic energy is

```math
\begin{align*}
    E &= \frac{1}{2V} ∫ |\bm{v}_{\text{BS}} + \bm{v}_{\text{f}}|^2 \, \mathrm{d}^3\bm{x}
    \\
    &= \frac{1}{2V} \left[
    Γ ∮_{\mathcal{C}} \left( \bm{ψ}_{\text{BS}} + 2 \bm{ψ}_{\text{f}} \right) ⋅ \mathrm{d}\bm{s}
    +
    ∫ |\bm{v}_{\text{f}}|^2 \, \mathrm{d}^3\bm{x}
    \right],
\end{align*}
```

where ``\bm{ψ}_{\text{BS}}`` is the streamfunction obtained from the Biot-Savart law.
Accordingly, when passing `external_streamfunction`, the streamfunction values attached to
vortex filaments (in `iter.ψs`) will actually contain
``\bm{ψ}_{\text{BS}} + 2 \bm{ψ}_{\text{f}}`` (note the factor 2) to enable the
computation of the full kinetic energy.
Moreover, calling [`Diagnostics.kinetic_energy_from_streamfunction`](@ref) with a
[`VortexFilamentSolver`](@ref) will automatically include an estimation of the kinetic
energy of the external velocity field (the volume integral above).

!!! note

    Whether this total kinetic energy is of any physical interest is a different question,
    so it can be reasonable to ignore the `external_streamfunction` parameter in order to
    obtain the energy associated to ``\bm{v}_{\text{BS}}`` only.
"""
function init(
        prob::VortexFilamentProblem{T}, scheme::TemporalScheme;
        alias_u0 = false,   # same as in OrdinaryDiffEq.jl
        dt::Real,
        dtmin::Real = 0.0,
        fast_term::FastBiotSavartTerm = LocalTerm(),
        refinement::RefinementCriterion = NoRefinement(),
        reconnect::ReconnectionCriterion = NoReconnections(),
        adaptivity::AdaptivityCriterion = NoAdaptivity(),
        callback::F = identity,
        external_velocity::ExtVel = nothing,
        external_streamfunction::ExtStf = nothing,
        timer = TimerOutput("VortexFilament"),
    ) where {
        F <: Function,
        ExtVel <: Union{Nothing, Function},
        ExtStf <: Union{Nothing, Function},
        T,
    }
    (; tspan,) = prob
    (; Ls,) = prob.p
    fs = convert(VectorOfVectors, prob.fs)
    vs_data = map(similar ∘ nodes, fs) :: AllFilamentVelocities
    vs = convert(VectorOfVectors, vs_data)
    ψs = similar(vs)
    fs_sol = alias_u0 ? fs : copy(fs)
    cache_reconnect = Reconnections.init_cache(reconnect, fs_sol, Ls)
    cache_bs = BiotSavart.init_cache(prob.p, fs_sol; timer)
    cache_timestepper = init_cache(scheme, fs, vs)

    # Wrap functions with the timer, so that timings are estimated each time the function is called.
    advect! = timer(advect_filaments!, "Advect filaments")
    rhs! = timer(update_values_at_nodes!, "Update values at nodes")
    callback_ = timer(callback, "Callback")

    if adaptivity !== NoAdaptivity() && !can_change_dt(scheme)
        throw(ArgumentError(lazy"temporal scheme $scheme doesn't support adaptibility; set `adaptivity = NoAdaptivity()` or choose a different scheme"))
    end

    time = TimeInfo(nstep = 0, t = first(tspan), dt = dt, dt_prev = dt)

    external_forcing = (
        velocity = external_velocity,
        streamfunction = external_streamfunction,
    )
    check_external_streamfunction(external_forcing, Ls)

    iter = VortexFilamentSolver(
        prob, fs_sol, vs, ψs, time, T(dtmin), refinement, adaptivity, cache_reconnect,
        cache_bs, cache_timestepper, fast_term, callback_, external_forcing,
        timer, advect!, rhs!,
    )

    finalise_step!(iter)

    iter
end

# Returns gradient of vector function f: ℝ³ ↦ ℝ³
function vector_gradient(f::F, x⃗::Vec3) where {F <: Function}
    ∇s = ntuple(Val(3)) do i
        # Gradient of i-th vector component.
        ForwardDiff.gradient(x -> f(x)[i], x⃗)
    end
    hcat(∇s...)  # gradient tensor (Aᵢⱼ = ∂ᵢfⱼ)
end

function curl(f::F, x⃗::Vec3) where {F <: Function}
    A = vector_gradient(f, x⃗)
    oftype(x⃗, (A[2, 3] - A[3, 2], A[3, 1] - A[1, 3], A[1, 2] - A[2, 1]))
end

# Check that the external streamfunction satisfies v = ∇ × ψ on a single arbitrary point.
# We check this using automatic differentiation of ψ at a given point.
function check_external_streamfunction(forcing::NamedTuple, Ls)
    (; velocity, streamfunction,) = forcing
    (velocity === nothing || streamfunction === nothing) && return nothing
    N = length(Ls)
    x⃗ = Vec3(ntuple(i -> 0.07 * i * Ls[i], Val(N)))  # arbitrary point in [0, L]³
    t = 0.0  # arbitrary time
    v⃗_actual = velocity(x⃗, t)
    v⃗_from_streamfunction = curl(y⃗ -> streamfunction(y⃗, t), x⃗)
    isapprox(v⃗_actual, v⃗_from_streamfunction) || @warn(
        """
        it seems like the external streamfunction doesn't match the external velocity.
        This may cause wrong energy estimates.
        """,
        v⃗_actual, v⃗_from_streamfunction,
    )
    nothing
end

function _update_values_at_nodes!(
        ::Val{:full},
        ::FastBiotSavartTerm,  # ignored in this case
        fields::NamedTuple{Names, NTuple{N, V}},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    ) where {Names, N, V <: VectorOfVectors}
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs)
end

# Compute slow component only.
# This is generally called in IMEX-RK substeps, where only the velocity (and not the
# streamfunction) is needed.
# We assume that the "slow" component is everything but LIA term when evolving the
# Biot-Savart law.
# This component will be treated explicitly by IMEX schemes.
function _update_values_at_nodes!(
        ::Val{:slow},
        ::LocalTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(false))
end

function _update_values_at_nodes!(
        ::Val{:slow},
        ::ShortRangeTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; shortrange = Val(false))
end

# Compute fast component only.
# This is generally called in IMEX-RK substeps, where only the velocity (and not the
# streamfunction) is needed.
# We assume that the "fast" component is the LIA term when evolving the Biot-Savart law.
# This component will be treated implicitly by IMEX schemes.
function _update_values_at_nodes!(
        ::Val{:fast},
        ::LocalTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; LIA = Val(:only))
end

function _update_values_at_nodes!(
        ::Val{:fast},
        ::ShortRangeTerm,
        fields::NamedTuple{(:velocity,)},
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    BiotSavart.compute_on_nodes!(fields, iter.cache_bs, fs; longrange = Val(false))
end

# This is the most general variant which should be called by timesteppers.
function update_values_at_nodes!(
        fields::NamedTuple, fs, t::Real, iter;
        component = Val(:full),  # compute slow + fast components by default
    )
    _update_values_at_nodes!(component, iter.fast_term, fields, fs, iter)
    _add_external_fields!(fields, iter.external_forcing, fs, t)
end

function _add_external_fields!(fields::NamedTuple, forcing::NamedTuple, fs, t)
    if haskey(fields, :velocity)
        _add_external_field!(fields.velocity, forcing.velocity, fs, t)
    end
    if haskey(fields, :streamfunction)
        # We multiply the external streamfunction by 2 to get the right kinetic energy.
        _add_external_field!(fields.streamfunction, forcing.streamfunction, fs, t; factor = 2)
    end
    fields
end

_add_external_field!(vs_all, ::Nothing, args...; kws...) = vs_all  # do nothing

function _add_external_field!(vs_all, vext::F, fs, time; factor = 1) where {F <: Function}
    @assert eachindex(vs_all) == eachindex(fs)
    for (f, vs) ∈ zip(fs, vs_all)
        for i ∈ eachindex(vs, f)
            @inbounds vs[i] = vs[i] + factor * vext(f[i], time)
        end
    end
    vs_all
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

function Reconnections.reconnect!(iter::VortexFilamentSolver)
    (; vs, ψs, fs, reconnect, to,) = iter
    fields = (vs, ψs)
    Reconnections.reconnect!(reconnect, fs; to) do f, i, mode
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
            # @info lazy"Removing filament with $(length(fs[i])) nodes"
            popat!(fs, i)
            popat!(vs, i)
            popat!(ψs, i)
        end
    end

    isempty(fs) && error("all vortices disappeared!")  # TODO nicer way to handle this?

    # Update velocities and streamfunctions to the next timestep (and first RK step).
    # Note that we only compute the streamfunction at full steps, and not in the middle of
    # RK substeps.
    # TODO make computation of ψ optional?
    fields = (velocity = vs, streamfunction = ψs,)

    # Note: here we always include the LIA terms, even when using IMEX or multirate schemes.
    # This must be taken into account by scheme implementations.
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

include("diagnostics.jl")

end
