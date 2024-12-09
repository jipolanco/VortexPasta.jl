export NoAdaptivity,
       AdaptBasedOnSegmentLength,
       AdaptBasedOnVelocity,
       MaximumTimestep

using ..BiotSavart: kelvin_wave_period

"""
    AdaptivityCriterion

Abstract type representing a temporal adaptivity criterion.

Implemented adaptivity criteria are:

- [`NoAdaptivity`](@ref): disables time adaptivity;

- [`AdaptBasedOnSegmentLength`](@ref): determines the timestep ``Δt`` based on the minimum
  distance ``ℓ_{\\min}`` between two filament nodes (``Δt ∝ ℓ_{\\min}^{-2}``);

- [`AdaptBasedOnVelocity`](@ref): determines the timestep ``Δt`` based on the maximum velocity
  ``v_{\\max}`` of filament nodes and on a predefined distance ``δ`` (``Δt = δ / v_{\\max}``).

## Combining multiple criteria

Adaptivity criteria can be combined using `|`. Example:

    adaptivity = AdaptBasedOnSegmentLength(1.4) | AdaptBasedOnVelocity(0.01)

As expected, the timestep ``Δt`` will be chosen so that it satisfies both criteria.
"""
abstract type AdaptivityCriterion end

"""
    maximum_displacement(::AdaptivityCriterion) -> Float64

Return the maximum node displacement ``δ_{\\text{crit}}`` allowed by the adaptivity criterion.

More precisely, ``δ_{\\text{crit}}`` controls the maximum displacement of a filament node
during a single timestep of duration ``Δt``.

This function is mainly relevant when using the [`AdaptBasedOnVelocity`](@ref) criterion.
For other criteria this returns `Inf`, meaning that they don't limit the maximum allowed
displacement.
"""
function maximum_displacement end

maximum_displacement(::AdaptivityCriterion) = Inf

"""
    NoAdaptivity <: AdaptivityCriterion
    NoAdaptivity()

Disable temporal adaptivity, leaving the timestep ``Δt`` constant.
"""
struct NoAdaptivity <: AdaptivityCriterion end

estimate_timestep(::NoAdaptivity, iter::AbstractSolver) = get_dt(iter)  # don't change current dt

@doc raw"""
    AdaptBasedOnSegmentLength <: AdaptivityCriterion
    AdaptBasedOnSegmentLength(γ::Float64)

Adapt timestep ``Δt`` based on the minimum distance ``δ`` between two filament nodes.

More precisely, the timestep is set to ``Δt = γ \, T_{\text{KW}}(δ)``, where ``γ`` is a
dimensionless factor to be chosen, and:

```math
T_{\text{KW}}(λ) = \frac{2 λ^2}{Γ} {\left[
    \ln\left( \frac{λ}{πa} \right) + \frac{1}{2} - (Δ + γ)
\right]}^{-1}
```

is the period of a Kelvin wave of wavelength ``λ``.
See [`ParamsBiotSavart`](@ref) for the definitions of the vortex parameters ``Γ``, ``a`` and
``Δ``.

This criterion is somewhat analogous to the CFL condition in grid-based
computations, and ``γ`` is the analogous of the maximum CFL number to be allowed.
As such, the right value of ``γ`` will depend on the chosen temporal scheme.

For example, the [`RK4`](@ref) scheme seems to require ``γ ≈ 1`` to remain stable.
"""
struct AdaptBasedOnSegmentLength <: AdaptivityCriterion
    γ :: Float64
end

function estimate_timestep(crit::AdaptBasedOnSegmentLength, iter::AbstractSolver)
    (; γ,) = crit
    (; prob, fs,)  = iter
    p = prob.p :: ParamsBiotSavart

    # 1. Determine minimum distance.
    # If we're using the chordal filament parametrisation, then the knot increments are
    # equal to the node distances (assuming knots were recently recomputed).
    parametrisation = Filaments.parametrisation(fs)
    δ = if parametrisation === Filaments.ChordalParametrisation()
        Filaments.minimum_knot_increment(fs)  # faster to evaluate
    else
        Filaments.minimum_node_distance(fs)
    end

    # 2. Determine timestep.
    dt = γ * kelvin_wave_period(p, δ)

    dt
end

@doc raw"""
    AdaptBasedOnVelocity <: AdaptivityCriterion
    AdaptBasedOnVelocity(δ_crit::Float64; safety_factor = 0.8)

Adapt timestep ``Δt`` based on the maximum velocity ``v_{\max}`` of filament nodes.

The objective is that the resulting maximum displacement ``δ_{\max} = v_{\max} Δt`` stays
below the given critical displacement `δ_crit`.

One application of this criterion is to ensure that reconnections happen in-between two
solver iterations (that is, to avoid that two filaments cross each other without
reconnecting). In this case, ``δ`` should be proportional to the chosen distance below
which reconnections are performed.

# Implementation details

This criterion is used in two ways:

1. A priori, to decide the ``Δt`` to be used in the next iteration given the velocities of
   filament nodes at the start of the iteration (at time ``t``).
   This ensures ``Δt ≤ δ_{\text{crit}} / v_{\max}`` where ``v_{\max}`` is computed at the
   start of the iteration.

2. A posteriori, after the actual node displacements ``δ`` from time ``t`` to time
   ``t + Δt`` have been computed (e.g. using some Runge--Kutta scheme).
   If some ``δ`` is larger than ``δ_{\text{crit}}``, then the iteration is recomputed after
   halving the timestep (``Δt → Δt/2``).
   In this case, the original displacements are thrown away (rejected).
   This process can be eventually repeated until the criterion is satisfied.

# Optional safety factor

The optional `safety_factor` should be `≤1`, which will further reduce the timestep
chosen a priori in step 1. This is to try to avoid too many rejected timesteps in step 2,
e.g. in case the actual advection velocity is larger than the velocity at the beginning of the
timestep.

# In combination with other criteria

In principle, using this criterion can lead to an infinite timestep when the velocities are
zero. For this reason, it's a good idea to combine this criterion with the
[`AdaptBasedOnSegmentLength`](@ref) or the [`MaximumTimestep`](@ref) criterion. For example:

    adaptivity = AdaptBasedOnVelocity(2.0) | AdaptBasedOnSegmentLength(0.9)
    adaptivity = AdaptBasedOnVelocity(2.0) | MaximumTimestep(0.01)

In fact, the second option is done automatically in [`init`](@ref) if only an
`AdaptBasedOnVelocity` is passed.
In that case, the maximum timestep is taken to be the `dt` passed to `init`.
"""
struct AdaptBasedOnVelocity <: AdaptivityCriterion
    δ :: Float64
    safety_factor :: Float64
    function AdaptBasedOnVelocity(δ; safety_factor = 0.8)
        safety_factor ≤ 1 || throw(ArgumentError(lazy"safety_factor should be ≤ 1; got safety_factor = $safety_factor"))
        new(δ, safety_factor)
    end
end

function Base.show(io::IO, crit::AdaptBasedOnVelocity)
    (; δ, safety_factor,) = crit
    print(io, lazy"AdaptBasedOnVelocity($δ; safety_factor = $safety_factor)")
end

maximum_displacement(crit::AdaptBasedOnVelocity) = crit.δ

function estimate_timestep(crit::AdaptBasedOnVelocity, iter::AbstractSolver)
    (; dt_prev,) = iter  # dt_prev is the latest used timestep
    T = typeof(dt_prev)
    v_max = maximum_vector_norm(iter.quantities.vL)
    dt_estimate::T = crit.safety_factor * crit.δ / v_max
    # Avoid increasing the dt too abruptly compared to the previous timestep.
    # This seems to help reduce the number of rejected timesteps, in case the dt was reduced
    # in the previous iteration due to the a posteriori displacements.
    min(dt_estimate, 2 * dt_prev)::T
end

"""
    MaximumTimestep <: AdaptivityCriterion
    MaximumTimestep(Δt_max::Float64)

Criterion ensuring that the timestep will be kept below a maximal value `Δt_max`.
"""
struct MaximumTimestep <: AdaptivityCriterion
    Δt :: Float64
end

estimate_timestep(crit::MaximumTimestep, iter::AbstractSolver) = crit.Δt

struct CombinedAdaptivityCriteria{
        Criteria <: Tuple{Vararg{AdaptivityCriterion}}
    } <: AdaptivityCriterion
    criteria :: Criteria
end

maximum_displacement(crit::CombinedAdaptivityCriteria) = minimum(maximum_displacement, crit.criteria)

# Pretty-printing
Base.show(io::IO, crit::CombinedAdaptivityCriteria) = _show_combined(io, crit.criteria...)

function _show_combined(io, c::AdaptivityCriterion, next...)
    print(io, c, " | ")
    _show_combined(io, next...)
end

_show_combined(io, c::AdaptivityCriterion) = show(io, c)

CombinedAdaptivityCriteria(args...) = CombinedAdaptivityCriteria(args)

Base.:|(a::AdaptivityCriterion, b::AdaptivityCriterion) =
    CombinedAdaptivityCriteria(a, b)
Base.:|(a::CombinedAdaptivityCriteria, b::AdaptivityCriterion) =
    CombinedAdaptivityCriteria(a.criteria..., b)

function estimate_timestep(crit::CombinedAdaptivityCriteria, iter::AbstractSolver)
    (; criteria,) = crit
    dts = map(c -> estimate_timestep(c, iter), criteria)
    @debug lazy"Estimated timesteps: $dts"
    min(dts...)  # choose the smallest estimated timestep
end

# If only the AdaptBasedOnVelocity criterion is passed to `init`, then combine it with a
# MaximumTimestep criterion. Otherwise, the dt can indefinitely increase when velocities are small.
possibly_add_max_timestep(crit::AdaptBasedOnVelocity, dt) = crit | MaximumTimestep(dt)
possibly_add_max_timestep(crit::AdaptivityCriterion, dt) = crit
