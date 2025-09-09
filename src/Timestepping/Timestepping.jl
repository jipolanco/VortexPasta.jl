"""
    Timestepping

Module defining timestepping solvers for vortex filament simulations.
"""
module Timestepping

export init, solve!, step!, VortexFilamentProblem, VortexFilamentSolver,
       save_checkpoint, load_checkpoint,
       inject_filament!,
       ShortRangeTerm, LocalTerm,
       MinimalEnergy,
       ParamsBiotSavart,                                          # from ..BiotSavart
       NoReconnections, ReconnectBasedOnDistance, ReconnectFast,  # from ..Reconnections
       reset_timer!  # from TimerOutputs

using ..VortexPasta: VortexPasta  # version()

using ..Containers: VectorOfVectors

using ..Filaments:
    Filaments,
    AbstractFilament,
    DiscretisationMethod,
    UnitTangent,
    CurvatureVector,
    Vec3,
    nodes,
    segments,
    knots,
    number_type,
    eltype_nested,
    close_filament!,
    filament_length,
    RefinementCriterion,
    NoRefinement

using ..FilamentIO: FilamentIO

using ..Reconnections:
    Reconnections,
    AbstractReconnectionCache,
    ReconnectionCriterion,
    NoReconnections,
    ReconnectBasedOnDistance,
    ReconnectFast,
    reconnect!

using ..BiotSavart:
    BiotSavart,
    ParamsBiotSavart,
    BiotSavartCache,
    VectorOfFilaments,
    VectorOfPositions,
    VectorOfVelocities,
    periods

using ..Forcing: Forcing, AbstractForcing, AbstractDissipation, NoForcing, NoDissipation,
    NormalFluidForcing, FourierBandForcing, FourierBandForcingBS, DissipationBS, SmallScaleDissipationBS

# Reuse same init, solve! and step! functions from the SciML ecosystem, to avoid clashes.
# See https://docs.sciml.ai/CommonSolve/stable/
import CommonSolve: init, solve!, step!
using LinearAlgebra: ×
using OhMyThreads: OhMyThreads, DynamicScheduler, SerialScheduler, tforeach

using Accessors: @delete  # allows to remove entry from (Named)Tuple
using TimerOutputs: TimerOutputs, TimerOutput, @timeit, reset_timer!

abstract type AbstractProblem end
abstract type AbstractSolver end

@enum SimulationStatus begin
    SUCCESS
    NO_VORTICES_LEFT
end

@enum SimulationMode begin
    MODE_DEFAULT
    MODE_MINIMAL_ENERGY
end

DefaultMode() = MODE_DEFAULT
MinimalEnergy() = MODE_MINIMAL_ENERGY

include("timesteppers/timesteppers.jl")
include("adaptivity.jl")
include("problem.jl")
include("simulation_stats.jl")

"""
    TimeInfo

Contains information on the current time and timestep of a solver.

Some useful fields are:

- `t::Float64`: current time;

- `dt::Float64`: timestep to be used in next iteration;

- `dt_prev::Float64` : timestep used in the last performed iteration;

- `nstep::Int`: number of timesteps performed until now;

- `nrejected::Int`: number of rejected iterations.

When using the [`AdaptBasedOnVelocity`](@ref) criterion, an iteration can be rejected if
the actual filament displacement is too large compared to what is imposed by the criterion.
In that case, the iteration will be recomputed with a smaller timestep (`dt → dt/2`).
"""
@kwdef mutable struct TimeInfo
    nstep     :: Int = 0
    nrejected :: Int = 0
    t       :: Float64 = 0
    dt      :: Float64 = 0
    dt_prev :: Float64 = 0
end

function Base.show(io::IO, time::TimeInfo)
    (; nstep, nrejected, t, dt, dt_prev,) = time
    print(io, "TimeInfo:")
    print(io, "\n - nstep     = ", nstep)
    print(io, "\n - t         = ", t)
    print(io, "\n - dt        = ", dt)
    print(io, "\n - dt_prev   = ", dt_prev)
    print(io, "\n - nrejected = ", nrejected)
end

function Base.summary(io::IO, time::TimeInfo)
    (; nstep, nrejected, t, dt, dt_prev,) = time
    print(io, "TimeInfo(nstep = $nstep, t = $t, dt = $dt, dt_prev = $dt_prev, nrejected = $nrejected)")
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

# Included fields

Some useful fields included in a `VortexFilamentSolver` are:

- `prob`: associated [`VortexFilamentProblem`](@ref) (including Biot–Savart parameters);

- `fs`: current state of vortex filaments in the system;

- `quantities`: a `NamedTuple` containing data on filament nodes at the current simulation
  timestep. See **Physical quantities** below for more details.

- `time`: a [`TimeInfo`](@ref) object containing information such as the current time and timestep;

- `stats`: a [`SimulationStats`](@ref) object containing information such as the total
  number of reconnections since the beginning of the simulation;

- `to`: a `TimerOutput`, which records the time spent on different functions;

- `cache_bs`: the Biot–Savart cache, which contains data from short- and
  long-range computations.

## Physical quantities

The `quantities` field is a `NamedTuple` containing instantaneous data on filament nodes.
For example, `quantities.vs` is the self-induced vortex velocity. To get the velocity of
filament `i` at node `j`, one can do `quantities.vs[i][j]`. Moreover, these quantities are
often interpolable, which means that one can do `quantities.vs[i](j, 0.5)` to get the
velocity in-between two nodes.

For convenience, if one has a `VortexFilamentSolver` `iter`, the shortcut `iter.vs` is
equivalent to `iter.quantities.vs` (this also applies to all other quantities).

### List of quantities

- `vs`: self-induced vortex velocity (due to Biot–Savart law). If enabled, this also
  includes the contribution of an `external_velocity` (see [`init`](@ref)).

- `ψs`: self-induced streamfunction on vortices. This is useful for estimating the kinetic
  energy associated to the induced velocity field. If enabled, this also includes the
  contribution of an `external_streamfunction` (see [`init`](@ref)).

If a normal fluid is present (e.g. by passing `forcing = NormalFluidForcing(...)` to [`init`](@ref)),
then the following quantities are also available:

- `vL`: actual filament velocity after mutual friction due to normal fluid (see [`NormalFluidForcing`](@ref)).

- `vn`: normal fluid velocity at vortex locations.

- `tangents`: unit tangents ``\\bm{s}'`` at vortex locations.

For convenience, if there is no normal fluid, then `vL` is defined as an alias of `vs`
(they're the same object).

If `forcing = FourierBandForcingBS(...)` is passed (see [`FourierBandForcingBS`](@ref)),
then the `vf` field contains the forcing term.
Similarly, if `dissipation` contains a dissipation term (e.g. [`DissipationBS`](@ref) or
[`SmallScaleDissipationBS`](@ref)), then the `vdiss` field contains the dissipation term.

"""
struct VortexFilamentSolver{
        T,
        Problem <: VortexFilamentProblem{T},
        Filaments <: VectorOfVectors{Vec3{T}, <:AbstractFilament{T}},
        Mode <: SimulationMode,
        Quantities <: NamedTuple,
        Refinement <: RefinementCriterion,
        Adaptivity <: AdaptivityCriterion,
        CacheReconnect <: AbstractReconnectionCache,
        CacheBS <: BiotSavartCache,
        CacheTimestepper <: TemporalSchemeCache,
        FastTerm <: FastBiotSavartTerm,
        Timer <: TimerOutput,
        Affect <: Function,
        AffectTime <: Function,
        Callback <: Function,
        ExternalFields <: NamedTuple,
        StretchingVelocity <: Union{Nothing, Function},
        Forcing <: AbstractForcing,
        ForcingCache,
        Dissipation <: AbstractDissipation,
        DissipationCache,
        AdvectFunction <: Function,
        RHSFunction <: Function,
    } <: AbstractSolver
    prob  :: Problem
    fs    :: Filaments
    mode  :: Mode
    quantities :: Quantities  # stored fields: (vs, ψs) + optional ones (vn, tangents, ...)
    time  :: TimeInfo
    stats :: SimulationStats{T}
    dtmin :: T
    refinement        :: Refinement
    adaptivity        :: Adaptivity
    reconnect         :: CacheReconnect
    cache_bs          :: CacheBS
    cache_timestepper :: CacheTimestepper
    fast_term         :: FastTerm
    LIA           :: Bool
    fold_periodic :: Bool
    affect!    :: Affect      # signature: affect!(iter)
    affect_t!  :: AffectTime  # signature: affect_t!(iter, t) where t is the current time
    callback :: Callback      # signature: callback(iter)
    external_fields    :: ExternalFields  # velocity and streamfunction forcing
    stretching_velocity :: StretchingVelocity
    forcing  :: Forcing  # forcing term added to the vortex velocity: vL = vs + vf (+ possibly other terms)
    forcing_cache  :: ForcingCache
    dissipation  :: Dissipation  # dissipation term added to the vortex velocity: vL = vs + vdiss (+ possibly other terms)
    dissipation_cache  :: DissipationCache
    to       :: Timer
    advect!  :: AdvectFunction  # function for advecting filaments with a known velocity
    rhs!     :: RHSFunction     # function for estimating filament velocities (and sometimes streamfunction) from their positions
end

function Base.show(io_in::IO, iter::VortexFilamentSolver)
    io = IOContext(io_in, :prefix => " │  ")  # this may be used by some `show` functions, for example for the forcing
    print(io, "VortexFilamentSolver with fields:")
    if iter.mode !== DefaultMode()
        print(io, "\n ├─ mode: ", iter.mode)
    end
    print(io, "\n ├─ prob: ", iter.prob)
    _print_summary(io, iter.fs; pre = "\n ├─ fs: ", post = "vortex filaments")
    _print_summary(io, iter.vL; pre = "\n ├─ vL: ", post = "vortex line velocity (vs + mutual friction + forcing)")
    _print_summary(io, iter.vs; pre = "\n ├─ vs: ", post = "self-induced superfluid velocity")
    if hasproperty(iter, :vn)
        _print_summary(io, iter.vn; pre = "\n ├─ vn: ", post = "normal fluid velocity")
    end
    if hasproperty(iter, :v_ns)
        _print_summary(io, iter.v_ns; pre = "\n ├─ v_ns: ", post = "slip velocity used in forcing")
    end
    if hasproperty(iter, :tangents)
        _print_summary(io, iter.tangents; pre = "\n ├─ tangents: ", post = "local unit tangent")
    end
    _print_summary(io, iter.ψs; pre = "\n ├─ ψs: ", post = "streamfunction vector")
    print(io, "\n ├─ time: ")
    summary(io, iter.time)
    print(io, "\n ├─ stats: ")
    summary(io, iter.stats)
    print(io, "\n ├─ dtmin: ", iter.dtmin)
    print(io, "\n ├─ refinement: ", iter.refinement)
    print(io, "\n ├─ adaptivity: ", iter.adaptivity)
    print(io, "\n ├─ reconnect: ", Reconnections.criterion(iter.reconnect))
    print(io, "\n ├─ fast_term: ", iter.fast_term)
    print(io, "\n ├─ LIA: ", iter.LIA)
    print(io, "\n ├─ fold_periodic: ", iter.fold_periodic)
    print(io, "\n ├─ cache_bs: ")
    summary(io, iter.cache_bs)
    print(io, "\n ├─ cache_timestepper: ")
    summary(io, iter.cache_timestepper)
    _maybe_print_function(io, "\n ├─ affect!:", iter.affect!)
    _maybe_print_function(io, "\n ├─ affect_t!:", iter.affect_t!)
    _maybe_print_function(io, "\n ├─ callback:", iter.callback)
    for (name, func) ∈ pairs(iter.external_fields)
        func === nothing || print(io, "\n ├─ external_fields.$name: Function ($func)")
    end
    if iter.stretching_velocity !== nothing
        print(io, "\n ├─ stretching_velocity: Function (", iter.stretching_velocity, ")")
    end
    if iter.forcing !== NoForcing()
        print(io, "\n ├─ forcing: ", iter.forcing)  # note: this draws a "subtree"
    end
    if iter.dissipation !== NoDissipation()
        print(io, "\n ├─ dissipation: ", iter.dissipation)  # note: this draws a "subtree"
    end
    # print(io, "\n ├─ advect!: Function")
    # print(io, "\n ├─ rhs!: Function")
    print(io, "\n └─ to: ")
    summary(io, iter.to)
end

# This is for convenience: allows doing e.g. `iter.t` instead of `iter.time.t` to get the
# current time.
@inline function Base.getproperty(iter::VortexFilamentSolver, name::Symbol)
    time = getfield(iter, :time)
    quantities = getfield(iter, :quantities)
    if hasproperty(time, name)
        getproperty(time, name)
    elseif haskey(quantities, name)
        quantities[name]
    else
        getfield(iter, name)
    end
end

function Base.propertynames(iter::VortexFilamentSolver, private::Bool = false)
    (fieldnames(typeof(iter))..., propertynames(iter.quantities, private)..., propertynames(iter.time, private)...)
end

# Default value of callback and affect arguments of init.
default_callback() = Returns(nothing)

# Wrap function with a TimerOutputs timer, returning an InstrumentedFunction which will
# report its running time each time it's called. We do _not_ do this if the function is the
# default callback function (which doesn't do anything), simply to avoid "noise" in the
# TimerOutputs outputs (e.g. calls to default affect! function when we haven't defined an
# actual function that does something interesting).
_maybe_wrap_function_with_timer(f::F, timer, name) where {F <: Function} = timer(f, name)::TimerOutputs.InstrumentedFunction
_maybe_wrap_function_with_timer(f::typeof(default_callback()), timer, name) = f  # don't wrap function if it doesn't do anything

# Unwrap InstrumentedFunction (for timings)
_printable_function(f::TimerOutputs.InstrumentedFunction) = f.func
_printable_function(f::Function) = f

function _maybe_print_function(io, prefix, f_to::Function)
    f = _printable_function(f_to)
    if f !== default_callback()
        m = methods(f)
        @assert length(m) ≥ 1  # usually there is only a single method
        print(io, prefix, " Function ", m[1])
    end
    nothing
end

get_dt(iter::VortexFilamentSolver) = iter.time.dt
get_t(iter::VortexFilamentSolver) = iter.time.t

@doc raw"""
    init(prob::VortexFilamentProblem, scheme::TemporalScheme; dt::Real, kws...) -> VortexFilamentSolver

Initialise vortex filament problem.

Returns a [`VortexFilamentSolver`](@ref) which can be advanced in time using
either [`step!`](@ref) or [`solve!`](@ref).

# Mandatory positional arguments

- `prob`: a [`VortexFilamentProblem`](@ref) containing the problem definition.

- `scheme`: a timestepping scheme (see [`TemporalScheme`](@ref)).

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

- `LIA = false`: if `true`, only use the local induction approximation (LIA) to advance
  vortex filaments, ignoring all non-local interactions. Note that reconnections may still
  be enabled.

- `fold_periodic = true`: if `true` (default), vortices will be recentred onto the main unit cell
  when using periodic boundary conditions. It may be convenient to disable this for
  visualisation purposes. This setting doesn't affect the results (velocity of each
  filament, reconnections, …), besides possible spatial translations of the filaments
  proportional to the domain period.

- `callback`: a function to be called at the end of each timestep. Its signature must be
  `callback(iter::VortexFilamentSolver)`. **Filaments should not be modified** by this
  function. In that case use `affect!` instead. See notes below for more details.

- `affect!`: similar to `callback`, but allows to modify filaments and other quantities at
  the end of each timestep. Its signature must be `affect!(iter::VortexFilamentSolver)`.
  See notes below for more details.

- `affect_t!`: similar to `affect!`, but allows to modify filaments and other quantities
  _within a timestep_ (e.g. within each Runge–Kutta substep). Its signature must be
  `affect_t!(iter::VortexFilamentSolver, time::Real)`. See notes below for more details.

- `timer = TimerOutput("VortexFilament")`: an optional `TimerOutput` for
  recording the time spent on different functions.

- `external_velocity` / `external_streamfunction`: allows to add an external velocity field to
  "force" the filaments. See "Adding an external velocity" below for more details.

- `stretching_velocity`: allows to add an external "stretching" velocity to the filaments.
  This can be considered as an external energy (or length) injection mechanism.
  However, it can also lead to spurious Kelvin waves.
  See "Adding a stretching velocity" below for more details.

- `forcing`: adds an extra term to the vortex velocities representing a forcing term which
  generally injects energy (but it can also dissipate). See the [`Forcing`](@ref) module for
  different options.

- `dissipation`: adds an extra term to the vortex velocities representing a dissipation term
  which generally dissipates energy. See the [`Forcing`](@ref) module for different options.

# Extended help

## Difference between `callback`, `affect!` and `affect_t!`

The difference between the `callback(iter)` and `affect!(iter)` functions is the time at
which they are called:

- the `affect!` function is called *before* performing Biot-Savart computations from the
  latest filament positions. In other words, the fields in `iter.quantities` are not
  synchronised with `iter.fs`, and it generally makes no sense to access them.
  Things like energy estimates will be incorrect if done in `affect!`. On the other hand,
  the `affect!` function **allows to modify `iter.fs`** before Biot–Savart computations are
  performed. In particular, to inject new filaments, call [`inject_filament!`](@ref) from
  within the `affect!` function.

- the `callback` function is called *after* performing Biot-Savart computations.
  This means that the fields in `iter.quantities` have been computed from the latest filament positions.
  However, one **must not modify `iter.fs`**, or otherwise the velocity `iter.vs` (which
  will be used at the next timestep) will no longer correspond to the filament positions.

In addition, one can also set an `affect_t!(iter, t)` function which can be used to modify
the state of the solver in a fine-grained manner, e.g. at each Runge–Kutta substep. Note
that this function must take an additional argument `t` (the current time) as a second argument.
This may be used, for instance, to apply a forcing which varies in time.

## Forcing / energy injection methods

There are multiple implemented methods for affecting and injecting energy onto the vortex system.
These methods can represent different physical mechanisms (or in fact be quite unphysical), and some
of them can make more sense than others when simulating e.g. zero-temperature superfluids.

### Adding an external velocity

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

#### Notes on energy estimation

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

### Adding a stretching velocity

One can set the **`stretching_velocity`** keyword argument to impose a stretching velocity
``\bm{v}_\text{L}(ξ, t)`` on filament locations ``\bm{s}(ξ, t)``.
This velocity will be parallel (and generally opposite) to the curvature vector
``\bm{s}'' = ρ \bm{N}``, where ``ρ`` is the curvature and ``\bm{N}`` the unit normal.

More precisely, the `stretching_velocity` parameter allows to set a filament velocity of the
form

```math
\bm{v}_{\text{L}}(ξ) = -v_{\text{L}}[ρ(ξ)] \, \bm{N}(ξ)
```

where ``v_{\text{L}}`` is a velocity magnitude which can optionally depend on the local
curvature value ``ρ``.
The `stretching_velocity` parameter should therefore be a function `v(ρ::Real) -> Real`.

Note that ``\bm{v} ⋅ \bm{s}''`` can be roughly interpreted as a local stretching rate
(this is actually true for the integrated quantity, see [`Diagnostics.stretching_rate`](@ref)).
Therefore, one may want ``v_{\text{L}}`` to be inversely proportional to the local
curvature ``ρ``, so that the local stretching rate is approximately independent of the
filament location ``ξ``.
To achieve this, one may set

```julia
stretching_velocity = ρ -> min(γ / ρ, v_max)
```

where `γ` is a constant (units ``T^{-1}``) setting the stretching magnitude, and `v_max` is
a maximum stretching velocity used to avoid very large velocities in low-curvature regions.
One could do something fancier and replace the `min` with a smooth regularisation such as

```julia
stretching_velocity = ρ -> -expm1(-ρ / ρ₀) * (γ / ρ)  # note: expm1(x) = exp(x) - 1
```

for some small curvature ``ρ₀``. The maximum allowed velocity will then be `vmax = γ / ρ₀`.

### Forcing via a normal fluid

Alternatively to the above approaches, one may impose a "normal" fluid velocity field to
modify, via a mutual friction force, the velocity of the vortex filaments. To do this, one
first needs to construct a [`NormalFluidForcing`](@ref) (see that link for definitions and examples)
representing a normal fluid velocity field. Then, the obtained object should be passed as
the optional `forcing` argument.

In principle, this method can be combined with the previous ones. The normal fluid forcing
will be applied _after_ all the other forcing methods.
In other words, the vortex velocity ``\bm{v}_{\text{s}}`` used to compute the mutual
friction velocity includes all the other contributions (i.e. external and stretching velocities).

There is also a [`FourierBandForcing`](@ref) type which behaves similarly to
`NormalFluidForcing` but is more localised in scale space, as it uses a band-pass filtered
version of the superfluid velocity induced by the vortices.

### Injecting filaments over time

Another way of injecting energy is simply by adding vortices to the simulation from time to
time. This can be achieved by using an `affect!` function. See [`inject_filament!`](@ref) for
some more details.

## Energy minimisation mode

It is also possible to run the simulation in "energy minimisation" mode, which can be useful
for finding a minimal energy state from a given initial condition. In this mode, filaments
will be advected not using the Biot–Savart velocity ``\bm{v}_{\text{s}}``, but using the velocity
``-\bm{s}' × \bm{v}_{\text{s}}`` where ``\bm{s}'`` is the local unit tangent vector.

This is because the functional derivative of the energy (per unit mass) with respect to the
vortex positions is:

```julia
\frac{δE}{δ\bm{s}} = Γ \bm{s}' × \bm{v}_{\text{s}}
```

In principle, the simulation should converge to a state that (locally) minimises the kinetic
energy. In some cases, e.g. when one starts with a single closed vortex, this is simply a
state where all vortices have disappeared.

Pass `mode = MinimalEnergy()` to enable this mode.
Forcing should be disabled to run in this mode.
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
        fold_periodic::Bool = true,
        LIA::Bool = false,
        callback::Callback = default_callback(),  # by default this is an empty function which just returns `nothing`
        affect!::Affect = default_callback(),
        affect_t!::AffectTime = default_callback(),
        external_velocity::ExtVel = nothing,
        external_streamfunction::ExtStf = nothing,
        stretching_velocity::StretchingVelocity = nothing,
        forcing::AbstractForcing = NoForcing(),
        dissipation::AbstractDissipation = NoDissipation(),
        mode::SimulationMode = DefaultMode(),
        timer = TimerOutput("VortexFilament"),
    ) where {
        Callback <: Function,
        Affect <: Function,
        AffectTime <: Function,
        ExtVel <: Union{Nothing, Function},
        ExtStf <: Union{Nothing, Function},
        StretchingVelocity <: Union{Nothing, Function},
        T,
    }
    (; tspan, state,) = prob
    (; Ls,) = prob.p
    fs = convert(VectorOfVectors, prob.fs)

    # Make sure the filaments are ready for interpolations and derivative estimations
    # for f in fs
    #     Filaments.update_coefficients!(f)
    # end

    vs = map(allocate_field_for_filament, fs)::VectorOfVectors
    ψs = similar(vs)
    tangents = map(similar ∘ nodes, fs)  # unit tangents (not interpolable)

    mode === MinimalEnergy() && forcing !== NoForcing() && throw(
        ArgumentError(
            "forcing should be disabled in MinimalEnergy mode"
        )
    )

    quantities_base = (; vs, ψs, tangents,)

    quantities_with_forcing = if forcing isa NormalFluidForcing
        (; quantities_base..., vL = similar(vs), vf = similar(vs), vn = similar(vs),)
    elseif forcing isa FourierBandForcing
        (; quantities_base..., v_ns = similar(vs), vf = similar(vs), vL = similar(vs),)  # we separately store vs and vL
    elseif forcing isa FourierBandForcingBS
        (; quantities_base..., vf = similar(vs), vL = similar(vs),)  # no normal fluid in this case
    elseif mode === MinimalEnergy() || dissipation !== NoDissipation()
        (; quantities_base..., vL = similar(vs),)  # vL and vs are different
    else
        (; quantities_base..., vL = vs,)  # vL is an alias to vs
    end

    quantities = if dissipation isa NoDissipation
        quantities_with_forcing
    else
        (; quantities_with_forcing..., vdiss = similar(vs),)  # include dissipation term vdiss
    end

    fs_sol = alias_u0 ? fs : copy(fs)

    cache_reconnect = Reconnections.init_cache(reconnect, fs_sol, Ls)
    cache_bs = BiotSavart.init_cache(prob.p, fs_sol; timer)
    cache_timestepper = init_cache(scheme, fs, vs)

    # Wrap functions with the timer, so that timings are estimated each time the function is called.
    advect! = timer(advect_filaments!, "Advect filaments")
    rhs! = timer(update_values_at_nodes!, "Update values at nodes")

    callback_ = _maybe_wrap_function_with_timer(callback, timer, "Callback")
    affect_ = _maybe_wrap_function_with_timer(affect!, timer, "Affect!")
    affect_t_ = _maybe_wrap_function_with_timer(affect_t!, timer, "Affect! (fine-grained)")

    adaptivity_ = possibly_add_max_timestep(adaptivity, dt)
    if adaptivity_ !== NoAdaptivity() && !can_change_dt(scheme)
        throw(ArgumentError(lazy"temporal scheme $scheme doesn't support adaptibility; set `adaptivity = NoAdaptivity()` or choose a different scheme"))
    end

    # It doesn't make much sense to combine the LIA and fast_term arguments, since there's
    # no RHS splitting to do when using LIA.
    # Therefore, if LIA is enabled, we only support the default fast_term = LocalTerm().
    if LIA && fast_term !== LocalTerm()
        throw(ArgumentError("currently, using LIA requires setting fast_term = LocalTerm()"))
    end

    if state === (;)  # empty state, new simulation
        time = TimeInfo(t = first(tspan), dt = dt, dt_prev = dt)
        stats = SimulationStats(T)
    else  # with state, restarting from checkpoint
        time = state.time
        if adaptivity === NoAdaptivity()
            time.dt = dt  # use constant timestep given as input
        end
        stats = state.stats::SimulationStats{T}
    end

    external_fields = (
        velocity = external_velocity,
        streamfunction = external_streamfunction,
    )
    check_external_streamfunction(external_fields, Ls)

    if stretching_velocity !== nothing
        # Check that the function accepts a real value of type T (local curvature) and
        # returns a real value (velocity magnitude).
        stretching_velocity(one(T))::Real
    end

    forcing_cache = Forcing.init_cache(forcing, cache_bs)
    dissipation_cache = Forcing.init_cache(dissipation, cache_bs)

    iter = VortexFilamentSolver(
        prob, fs_sol, mode, quantities, time, stats, T(dtmin), refinement, adaptivity_, cache_reconnect,
        cache_bs, cache_timestepper, fast_term, LIA, fold_periodic, affect_, affect_t_, callback_, external_fields,
        stretching_velocity, forcing, forcing_cache, dissipation, dissipation_cache,
        timer, advect!, rhs!,
    )

    # Verify signature of callback and affect functions
    applicable(callback, iter) || throw(ArgumentError("`callback` function should be callable as `f(iter::VortexFilamentSolver)`"))
    applicable(affect!, iter) || throw(ArgumentError("`affect!` function should be callable as `f(iter::VortexFilamentSolver`)"))
    applicable(affect_t!, iter, time.t) || throw(ArgumentError("`affect_t!` function should be callable as `f(iter::VortexFilamentSolver, t::Real)`"))

    status = finalise_step!(iter)
    status == SUCCESS || error("reached status = $status at initialisation")

    iter
end

function allocate_field_for_filament(f::AbstractFilament)
    allocate_field_for_filament(typeof(f), f)  # by default, return a Filament
end

# This is to allocate a vector to hold quantities such as velocity, ...
# We currently represent them as filaments because we want them to be interpolable,
# as this is convenient for some diagnostics.
# For now we don't require any derivatives.
function allocate_field_for_filament(::Type{<:AbstractFilament}, f::AbstractFilament)
    Filaments.similar_filament(f; nderivs = Val(0), offset = zero(eltype(f)))
end

# Variant that returns a PaddedVector instead of a Filament
function allocate_field_for_filament(::Type{<:AbstractVector}, f::AbstractFilament)
    similar(Filaments.nodes(f))
end

function fields_to_resize(iter::VortexFilamentSolver)
    (; quantities,) = iter
    if quantities.vs === quantities.vL  # they're aliased
        values(@delete quantities.vL)  # resize all fields except vL (aliased with vs)
    else
        values(quantities)  # resize all fields
    end
end

function fields_to_interpolate(iter::VortexFilamentSolver)
    (; quantities,) = iter
    qs = @delete quantities.tangents  # don't interpolate tangents
    if qs.vs === qs.vL  # they're aliased
        values(@delete qs.vL)  # interpolate all fields except vL (aliased with vs)
    else
        values(qs)  # interpolate all fields (not sure we need all of them...)
    end
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
    v⃗_from_streamfunction = Filaments.curl(y⃗ -> streamfunction(y⃗, t), x⃗)
    isapprox(v⃗_actual, v⃗_from_streamfunction) || @warn(
        """
        it seems like the external streamfunction doesn't match the external velocity.
        This may cause wrong energy estimates.
        """,
        v⃗_actual, v⃗_from_streamfunction,
    )
    nothing
end

function refine!(f::AbstractFilament, refinement::RefinementCriterion)
    nref = Filaments.refine!(f, refinement)  # we assume update_coefficients! was already called
    n = 0
    nmax = 1  # maximum number of extra refinement iterations (to avoid infinite loop...)
    while nref !== (0, 0) && n < nmax
        n += 1
        @debug "Added/removed $nref nodes"
        if !Filaments.check_nodes(Bool, f)
            # This may happen if the refinement removed nodes, and now the number of nodes
            # is too small (usually < 3) to represent the filament.
            # It doesn't make sense to keep refining the filament (and it will fail
            # anyways).
            # The filament will be removed in finalise_step!.
            break
        end
        nref = Filaments.refine!(f, refinement)
    end
    n
end

function after_advection!(iter::VortexFilamentSolver)
    (; fs, prob, refinement, fold_periodic, stats, to,) = iter
    # Perform reconnections, possibly changing the number of filaments.
    rec = reconnect!(iter)
    @debug lazy"Number of reconnections: $(rec.reconnection_count)"
    stats.reconnection_count += rec.reconnection_count
    stats.reconnection_length_loss += rec.reconnection_length_loss
    stats.filaments_removed_count += rec.filaments_removed_count
    stats.filaments_removed_length += rec.filaments_removed_length
    stats.reconnection_passes += rec.npasses
    L_fold = periods(prob.p)  # box size (periodicity)
    fields = fields_to_resize(iter)
    after_advection!(fs, fields; L_fold, refinement, fold_periodic,)
end

# This should be called right after positions have been updated by the timestepper.
# Here `fields` is a tuple containing the fields associated to each filament node.
# Usually this is `(vs, ψs)`, which contain the velocities and streamfunction values at all
# filament nodes.
function after_advection!(fs, fields::Tuple; kws...)
    @inbounds for i ∈ eachindex(fs)
        fields_i = map(vs -> @inbounds(vs[i]), fields)  # typically this is (vs[i], ψs[i])
        _after_advection!(fs[i], fields_i; kws...)
    end
    fs
end

function _after_advection!(f::AbstractFilament, fields::Tuple; L_fold, refinement, fold_periodic,)
    if fold_periodic && Filaments.fold_periodic!(f, L_fold)
        Filaments.update_coefficients!(f)  # only called if nodes were modified by fold_periodic!
    end
    refinement_steps = refine!(f, refinement)
    if refinement_steps > 0  # refinement was performed
        map(vs -> resize!(vs, length(f)), fields)
    end
    f
end

"""
    solve!(iter::VortexFilamentSolver)

Advance vortex filament solver to the ending time.

See also [`step!`](@ref) for advancing one step at a time.
"""
function solve!(iter::VortexFilamentSolver)
    (; time, adaptivity,) = iter
    t_end = iter.prob.tspan[2]
    while time.t < t_end
        if adaptivity === NoAdaptivity()  # constant timestep
            if isapprox(time.t, t_end; atol = time.dt / 1000)
                break  # stop simulation if we're very close to the end time
            end
        elseif can_change_dt(iter.cache_timestepper)
            # Try to finish exactly at t = t_end.
            # Note: we don't do this when using a constant timestep.
            time.dt = min(time.dt, t_end - time.t)
        end
        status = step!(iter)
        if status != SUCCESS
            @warn lazy"reached status = $status. Stopping the simulation."
            break
        end
    end
    iter
end

# Returns the maximum |v⃗| from a VectorOfVectors.
# We mainly use it to get the maximum velocity norm among all filament nodes.
function maximum_vector_norm(vs::VectorOfVectors{<:Vec3})
    T = number_type(vs)
    @assert T <: AbstractFloat
    v²_max = zero(T)
    for vnodes ∈ vs, v⃗ ∈ vnodes
        v²_max = max(v²_max, sum(abs2, v⃗)::T)
    end
    sqrt(v²_max)
end

"""
    step!(iter::VortexFilamentSolver) -> SimulationStatus

Advance solver by a single timestep.

Returns a `SimulationStatus`. Currently, this can be:

- `SUCCESS`,
- `NO_VORTICES_LEFT`: all vortices have been removed from the system; the simulation should
  be stopped.
"""
function step!(iter::VortexFilamentSolver{T}) where {T}
    (; fs, quantities, time, adaptivity, dtmin, advect!, rhs!,) = iter
    (; ψs, vL,) = quantities
    vL_start = ψs  # reuse ψs to store velocity at start of timestep (needed by RK schemes)
    copyto!(vL_start, vL)
    t_end = iter.prob.tspan[2]
    if time.dt < dtmin && time.t + time.dt < t_end
        error(lazy"current timestep is too small ($(time.dt) < $(dtmin)). Stopping.")
    end

    # Maximum allowed displacement in a single timestep (can be Inf).
    δ_crit::T = maximum_displacement(adaptivity)

    while true
        # Note: the timesteppers assume that vL already contains the filament velocities
        # computed at the start of the current timestep.
        update_velocities!(vL, rhs!, advect!, iter.cache_timestepper, iter)  # timestepping

        # Compute actual maximum displacement.
        v_max = maximum_vector_norm(vL)::T
        δ_max = v_max * T(time.dt)

        if δ_max > δ_crit
            # Reject velocities and reduce timestep by half.
            time.nrejected += 1
            time.dt /= 2
            copyto!(vL, vL_start)  # reset velocities to the beginning of timestep
        else
            break
        end
    end

    advect!(fs, vL, time.dt)
    after_advection!(iter)
    time.t += time.dt
    time.nstep += 1
    status = finalise_step!(iter)

    status
end

# Here fields is usually (vs, ψs, vL, ...)
@inline function reconnect_callback((fs, fields), f, i, mode::Symbol)
    if mode === :removed
        @debug lazy"Filament was removed at index $i"
        # @assert i ≤ lastindex(fields[1]) == lastindex(fs) + 1
        map(vs -> popat!(vs, i), fields)
    elseif mode === :appended
        @debug lazy"Filament was appended at index $i"
        @assert f === fs[i]
        # @assert i == lastindex(fields[1]) + 1
        map(fields) do vs
            vf = allocate_field_for_filament(eltype(vs), f)
            push!(vs, vf)
        end
    elseif mode === :modified
        @debug lazy"Filament was modified at index $i"
        @assert f === fs[i]
        # @assert i ≤ lastindex(fs) == lastindex(fields[1])
        map(vs -> resize!(vs[i], length(f)), fields)
    end
    nothing
end

function Reconnections.reconnect!(iter::VortexFilamentSolver)
    # Here vL is the velocity that was used to perform the latest filament displacements (in advect!).
    (; vL, fs, reconnect, to,) = iter
    fields = fields_to_resize(iter)
    Reconnections.reconnect!(reconnect, fs, vL; to) do f, i, mode
        reconnect_callback((fs, fields), f, i, mode)
    end
end

# Called whenever filament positions have just been initialised or updated.
# Returns a SimulationStatus.
function finalise_step!(iter::VortexFilamentSolver)
    (; fs, time, stats, adaptivity, rhs!,) = iter
    (; vs, vL, ψs,) = iter.quantities

    @assert eachindex(fs) == eachindex(vL) == eachindex(vs)

    # Check if there are filaments to be removed (typically with ≤ 3 discretisation points,
    # but this depends on the actual discretisation method). Note that filaments may also
    # be removed during reconnections for the same reason.
    for i ∈ reverse(eachindex(fs))
        f = fs[i]
        if !Filaments.check_nodes(Bool, f)
            # Remove filament and its associated vectors of velocities/streamfunctions.
            @debug lazy"Removing filament with $(length(f)) nodes"
            # Compute vortex length using the straight segment approximation (no quadratures),
            # since the filament doesn't accept a continuous description at this point.
            Filaments.close_filament!(f)  # needed before calling filament_length
            local Lfil = filament_length(f; quad = nothing)
            stats.filaments_removed_count += 1
            stats.filaments_removed_length += Lfil
            popat!(fs, i)
            for field ∈ fields_to_resize(iter)
                popat!(field, i)
            end
        end
    end

    isempty(fs) && return NO_VORTICES_LEFT

    iter.affect!(iter)

    # Update velocities and streamfunctions to the next timestep (and first RK step).
    # Note that we only compute the streamfunction at full steps, and not in the middle of
    # RK substeps.
    # TODO: make computation of ψ optional?
    fields = (velocity = vL, streamfunction = ψs,)

    # Note: here we always include the LIA terms, even when using IMEX or multirate schemes.
    # This must be taken into account by scheme implementations.
    rhs!(fields, fs, time.t, iter; component = Val(:full))

    # Update coefficients in case we need to interpolate fields, e.g. for diagnostics or
    # reconnections. We make sure we use the same parametrisation (knots) of the filaments
    # themselves.
    # TODO: perform reconnections after this?
    for i ∈ eachindex(fs)
        local ts = Filaments.knots(fs[i])
        foreach(fields_to_interpolate(iter)) do qs
            Filaments.update_coefficients!(qs[i]; knots = ts)
        end
    end

    time.dt_prev = time.dt
    time.dt = estimate_timestep(adaptivity, iter)  # estimate dt for next timestep
    iter.callback(iter)

    SUCCESS
end

"""
    inject_filament!(iter::VortexFilamentSolver, f::AbstractFilament)

Inject a filament onto a simulation.

This function **must** be called from within the `affect!` function attached to the solver
(see [`init`](@ref) for details). Otherwise, the solver will likely use garbage values for
the velocity of the injected filament, leading to wrong results or worse.

# Basic usage

A very simplified (and incomplete) example of how to inject filaments during a simulation:

```julia
# Define function that will be called after each simulation timestep
function my_affect!(iter::VortexFilamentSolver)
    f = Filaments.init(...)  # create new filament
    inject_filament!(iter, f)
    return nothing
end

# Initialise and run simulation
prob = VortexFilamentProblem(...)
iter = init(prob, RK4(); affect! = my_affect!, etc...)
solve!(iter)
```
"""
function inject_filament!(iter::VortexFilamentSolver, f::AbstractFilament)
    (; fs,) = iter
    # 1. Add filament to list of filaments
    push!(fs, f)
    # 2. Add space to store values on filaments (velocities, etc.)
    for us ∈ fields_to_resize(iter)
        @assert !isempty(us)
        u_new = similar(first(us), length(f))
        push!(us, u_new)
    end
    nothing
end

include("advect.jl")  # advect_filaments! (a.k.a. iter.advect!)
include("rhs.jl")     # update_values_at_nodes! (a.k.a. iter.rhs!)
include("forcing.jl")
include("diagnostics.jl")
include("checkpoint.jl")

end

