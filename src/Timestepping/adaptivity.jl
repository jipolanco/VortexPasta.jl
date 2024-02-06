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

"""
    AdaptBasedOnVelocity <: AdaptivityCriterion
    AdaptBasedOnVelocity(δ::Float64)

Adapt timestep ``Δt`` based on the maximum velocity ``v_{\\max}`` of filament nodes and on
the given distance ``δ``.

The timestep is set to ``Δt = δ / v_{\\max}``.

Note that, in principle, using this criterion can lead to an infinite timestep when the
velocities are zero. For this reason, it's a good idea to combine this criterion with the
[`MaximumTimestep`](@ref) criterion. For example:

    adaptivity = AdaptBasedOnVelocity(2.0) | MaximumTimestep(0.01)

In fact, this is done automatically in [`init`](@ref) if only an `AdaptBasedOnVelocity` is passed.
In that case, the maximum timestep is taken to be the `dt` passed to `init`.

One application of this criterion is to ensure that reconnections happen in-between two
solver iterations (that is, to avoid that two filaments cross each other without
reconnecting). In this case, ``δ`` should be proportional to the chosen distance below
which reconnections are performed.
"""
struct AdaptBasedOnVelocity <: AdaptivityCriterion
    δ :: Float64
end

function estimate_timestep(crit::AdaptBasedOnVelocity, iter::AbstractSolver)
    (; δ,) = crit
    (; vs,)  = iter
    T = eltype(eltype(eltype(vs)))
    @assert T <: AbstractFloat
    v²_max = maximum(vs; init = zero(T)) do vnodes
        maximum(v⃗ -> sum(abs2, v⃗), vnodes)  # maximum squared velocity norm among the nodes of a single filament
    end
    v_max = sqrt(v²_max)
    δ / v_max
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
