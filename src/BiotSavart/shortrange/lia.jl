using ..Filaments: CurvatureBinormal, UnitTangent
using LinearAlgebra: normalize
using StaticArrays: SVector

"""
    local_self_induced_velocity(
        f::AbstractFilament, i::Int, [prefactor::Real];
        a::Real, [Γ::Real],
        Δ = 0.25, quad = nothing, fit_circle = false,
        segment_fraction = nothing,
    )

Compute local self-induced velocity of filament node `f[i]`.

This corresponds to the LIA term (localised induction approximation).

This is the same as `local_self_induced(Velocity(), f, i, [prefactor]; kwargs...)`.
See also [`local_self_induced`](@ref).
"""
local_self_induced_velocity(args...; kws...) =
    local_self_induced(Velocity(), args...; kws...)

# Prefactor was not given as an argument, compute it.
local_self_induced(q::OutputField, f, i; Γ, kws...) =
    local_self_induced(q, f, i, Γ / 4π; kws...)

"""
    local_self_induced(
        q::OutputField, f::AbstractFilament, i::Int, [prefactor::Real];
        a::Real, Δ::Real = 0.25, segment_fraction = nothing, quad = nothing, [Γ],
    )

Compute localised induction approximation (LIA) term at node `f[i]`.

Possible output fields to be passed as first argument are `Velocity()` and `Streamfunction()`.

One must pass either the circulation `Γ` as a keyword argument, or a precomputed
`prefactor` which is usually equal to `Γ / 4π` (both for velocity and streamfunction).

# Mandatory arguments

- the vortex core size `a` as a keyword argument;

- either the vortex circulation `Γ` as a keyword argument, or the precomputed
  prefactor `Γ / 4π` as a positional argument.

# Optional arguments

- the optional parameter `Δ` sets the effect of the vorticity profile (see
  [`ParamsBiotSavart`](@ref) for details);

- the optional parameter `segment_fraction` can be used to indicate that the LIA term should
  only be computed over a fraction of the segments neighbouring the node `f[i]`. In that case,
  it should be a real number in ``(0, 1]`` (where 1 corresponds to integrating over the full
  segments, which the default);

- a quadrature rule can be passed via `quad`, which can improve the estimation
  of the LIA term and the stability of the solver (even a 1-point quadrature
  rule can importantly improve stability!);

- if `quad = nothing`, one may set `fit_circle = true` to estimate the binormal
  vector by fitting a circle passing through 3 neighbouring nodes (as done in
  Schwarz PRB 1985), instead of using local derivatives. For now, this is only implemented
  for `q = Velocity()`, and it may be removed in the future.

"""
function local_self_induced(
        q::OutputField, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real = 0.25, quad = nothing, kws...,
    )
    _local_self_induced(q, to_quadrature(quad), f, i, prefactor; a, Δ, kws...)
end

to_quadrature(quad::AbstractQuadrature) = quad
to_quadrature(::Nothing) = NoQuadrature()

function _local_self_induced(
        ::Velocity, ::NoQuadrature, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real, fit_circle = false, segment_fraction = nothing, kws...,
    )
    ℓ₋ = norm(f[i] - f[i - 1])
    ℓ₊ = norm(f[i + 1] - f[i])
    γ = something(segment_fraction, true)  # true in the sense of 1
    β = prefactor * (log(2 * γ * sqrt(ℓ₋ * ℓ₊) / a) - Δ)
    if fit_circle
        # Fit circle passing through the 3 points.
        # See https://en.wikipedia.org/wiki/Circumscribed_circle#Higher_dimensions
        A = @inbounds f[i - 1]
        B = @inbounds f[i]
        C = @inbounds f[i + 1]
        a = A - B
        b = C - B
        a2 = sum(abs2, a)
        b2 = sum(abs2, b)
        ab2 = sum(abs2, a - b)
        (2 * β / sqrt(a2 * b2 * ab2)) * (b × a)
    else
        # Evaluate derivatives at node `i`.
        Ẋ = f[i, Derivative(1)]  # ∂f/∂t at node i
        Ẍ = f[i, Derivative(2)]
        # Note that the derivatives are wrt the parameter `t` and not exactly wrt
        # the arc length `ξ`, hence the extra `norm(Ẋ)^3` factor wrt the literature.
        # Usually, `norm(Ẋ)` should be actually quite close to 1 for chordal
        # parametrisation.
        β / norm(Ẋ)^3 * (Ẋ × Ẍ)
    end
end

lia_integration_limits(::Nothing) = (nothing, nothing)
lia_integration_limits(γ::Real) = ((one(γ) - γ, one(γ)), (zero(γ), γ))

nonlia_integration_limits(::Nothing) = (nothing, nothing)
nonlia_integration_limits(γ::Real) = ((zero(γ), one(γ) - γ), (γ, one(γ)))

# Alternative estimation using quadratures.
@inline function _local_self_induced(
        ::Velocity, quad::AbstractQuadrature, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real, segment_fraction::Union{Nothing, Real} = nothing,
        kws...,
    )
    lims = lia_integration_limits(segment_fraction)
    ℓ₋² = integrate(f, i - 1, quad; limits = lims[1]) do f, j, ζ
        sum(abs2, f(j, ζ, Derivative(1)))
    end
    ℓ₊² = integrate(f, i, quad; limits = lims[2]) do f, j, ζ
        sum(abs2, f(j, ζ, Derivative(1)))
    end
    b⃗ = f[i, CurvatureBinormal()]
    ℓ = sqrt(sqrt(ℓ₋² * ℓ₊²))
    β = prefactor * (log(2 * ℓ / a) - Δ)
    β * b⃗
end

function _local_self_induced(
        ::Streamfunction, ::NoQuadrature,
        f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real, segment_fraction = nothing, kws...,
    )
    ℓ₋ = norm(f[i] - f[i - 1])
    ℓ₊ = norm(f[i + 1] - f[i])
    γ = something(segment_fraction, true)  # true in the sense of 1
    t̂ = f[i, UnitTangent()]
    β = 2 * prefactor * (log(2 * γ * sqrt(ℓ₋ * ℓ₊) / a) + 1 - Δ)  # note: prefactor = Γ/4π (hence the 2 in front)
    β * t̂
end

@inline function _local_self_induced(
        ::Streamfunction, quad::AbstractQuadrature,
        f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real,
        segment_fraction::Union{Nothing, Real} = nothing,
        kws...,
    )
    lims = lia_integration_limits(segment_fraction)
    ℓ₋ = integrate(f, i - 1, quad; limits = lims[1]) do f, j, ζ
        norm(f(j, ζ, Derivative(1)))
    end
    ℓ₊ = integrate(f, i, quad; limits = lims[2]) do f, j, ζ
        norm(f(j, ζ, Derivative(1)))
    end
    # Use local (non-averaged) tangent at point of interest.
    # We want to make sure that the result is exactly tangent to the curve at the node.
    t̂ = f[i, UnitTangent()]
    # Note: the +1 coefficient is required for energy conservation.
    # It is required so that the resulting energy follows Hamilton's equation (at least for
    # the case of a vortex ring), and it has been verified in many cases that it improves
    # the effective energy conservation.
    β = 2 * prefactor * (log(2 * sqrt(ℓ₋ * ℓ₊) / a) + 1 - Δ)  # note: prefactor = Γ/4π (hence the 2 in front)
    β * t̂
end
