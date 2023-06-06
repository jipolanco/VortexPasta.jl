export NoAdaptivity,
       AdaptBasedOnSegmentLength,
       AdaptBasedOnVelocity

using LinearAlgebra: norm

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

"""
    AdaptBasedOnSegmentLength <: AdaptivityCriterion
    AdaptBasedOnSegmentLength(γ::Float64)

Adapt timestep ``Δt`` based on the minimum distance ``ℓ_{\\min}`` between two filament nodes.

More precisely, the timestep is set to ``Δt = γ ℓ_{\\min}² / β(ℓ_{\\min})``, where ``γ`` is a
dimensionless factor to be chosen. Here ``β(ℓ) = Γ / 4π \\left[ \\ln(2ℓ_{\\min} / a) - Δ
\\right]`` (see [`ParamsBiotSavart`](@ref) for the definitions of ``Γ``, ``a``
and ``Δ``).

Note that this criterion is based on simple dimensional analysis, assuming that
``Δt`` must only depend on the minimum segment length ``ℓ_{\\min}`` and on the vortex
circulation ``Γ`` (or rather on ``β``, which is proportional to ``Γ`` and also
includes the effect of the vortex core size and profile).

This criterion is somewhat analogous to the CFL condition in grid-based
computations, and ``γ`` is the analogous of the maximum CFL number to be allowed.
As such, the right value of ``γ`` will depend on the chosen temporal scheme.

Empirically, the following values of ``γ`` seem to provide stability and good precision:

- ``γ = 1.8`` for [`RK4`](@ref),
- ``γ = 1.2`` for [`SSPRK33`](@ref),
- ``γ = 0.12`` for [`Euler`](@ref).

This basically shows that among these methods, `RK4` should be preferred as it
allows for a larger timestep for a given number of function evaluations. On the
other hand, `Euler` should be always avoided!
"""
struct AdaptBasedOnSegmentLength <: AdaptivityCriterion
    γ :: Float64
end

function estimate_timestep(crit::AdaptBasedOnSegmentLength, iter::AbstractSolver)
    (; γ,) = crit
    (; prob, fs,)  = iter
    p = prob.p :: ParamsBiotSavart

    # 1. Determine minimum distance ℓ_min between filament nodes.
    # TODO can we avoid iterating multiple times over the filaments?
    T = eltype(p)
    @assert T <: AbstractFloat
    ℓ_min = convert(T, Inf)
    @inbounds for f ∈ fs
        ts = knots(f)
        for i ∈ eachindex(segments(f))
            ℓ::T = ts[i + 1] - ts[i]  # assume that knots roughly correspond to segment lengths
            # @assert ℓ > 0
            ℓ_min = min(ℓ_min, ℓ)
        end
    end

    # 2. Determine timestep.
    (; Γ, a, Δ,) = p.common
    β = (Γ / 4π) * (log(2 * ℓ_min / a) - Δ)
    dt = γ * ℓ_min^2 / β

    dt
end

"""
    AdaptBasedOnVelocity <: AdaptivityCriterion
    AdaptBasedOnVelocity(δ::Float64)

Adapt timestep ``Δt`` based on the maximum velocity ``v_{\\max}`` of filament nodes and on
the given distance ``δ`.

The timestep is set to ``Δt = δ / v_{\\max}``.

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
    vmax = maximum(vs) do vnodes
        maximum(norm, vnodes)  # maximum velocity norm among the nodes of a single filament
    end
    δ / vmax
end

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
    minimum(c -> estimate_timestep(c, iter), criteria)  # choose the smallest estimated timestep
end
