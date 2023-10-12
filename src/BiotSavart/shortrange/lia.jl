using ..Filaments: CurvatureBinormal, UnitTangent
using LinearAlgebra: normalize
using StaticArrays: SVector

"""
    local_self_induced_velocity(
        f::AbstractFilament, i::Int, [prefactor::Real];
        a::Real, [Γ::Real],
        Δ = 0.25, quad = nothing, fit_circle = false,
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
        a::Real, Δ::Real = 0.25, quad = nothing, [Γ],
    )

Compute localised induction approximation (LIA) term at node `f[i]`.

Possible output fields to be passed as first argument are `Velocity()` and `Streamfunction()`.

One must pass either the circulation `Γ` as a keyword argument, or a precomputed
`prefactor` which is usually equal to `Γ / 4π` (both for velocity and streamfunction).

## Mandatory arguments

- the vortex core size `a` as a keyword argument;

- either the vortex circulation `Γ` as a keyword argument, or the precomputed
  prefactor `Γ / 4π` as a positional argument.

## Optional arguments

- the optional parameter `Δ` sets the effect of the vorticity profile (see
  [`ParamsBiotSavart`](@ref) for details);

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
    _local_self_induced(q, quad, f, i, prefactor; a, Δ, kws...)
end

function _local_self_induced(
        ::Velocity, quad::Nothing, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real, fit_circle = false,
    )
    ℓ₋ = norm(f[i] - f[i - 1])
    ℓ₊ = norm(f[i + 1] - f[i])
    β = prefactor * (log(2 * sqrt(ℓ₋ * ℓ₊) / a) - Δ)
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
        # Usually, `norm(Ẋ)` should be actually quite close to 1 since the
        # parametrisation is a rough approximation of the arc length.
        β / norm(Ẋ)^3 * (Ẋ × Ẍ)
    end
end

# Alternative estimation using quadratures.
# It seems to improve accuracy and stability (tested with vortex ring example and Kelvin waves).
function _local_self_induced(
        ::Velocity, quad::AbstractQuadrature, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real,
        fit_circle = false,  # ignored
    )
    ts = knots(f)
    ℓ₋ = integrate(f, i - 1, quad) do f, j, ζ
        norm(f(j, ζ, Derivative(1)))
    end
    ℓ₊ = integrate(f, i, quad) do f, j, ζ
        norm(f(j, ζ, Derivative(1)))
    end
    # Estimate the scaled binormal vector b⃗ = ρ b̂, where ρ is the curvature and b̂ = t̂ × n̂.
    b⃗₋ = integrate(f, i - 1, quad) do f, j, ζ
        f(j, ζ, CurvatureBinormal())
    end
    b⃗₊ = integrate(f, i, quad) do f, j, ζ
        f(j, ζ, CurvatureBinormal())
    end
    b⃗ = (b⃗₋ + b⃗₊) ./ (ts[i + 1] - ts[i - 1])  # average on [i - 1, i + 1]
    β = prefactor * (log(2 * sqrt(ℓ₋ * ℓ₊) / a) - Δ)
    β * b⃗
end

# TODO
# - Implement variant with no quadratures?
function _local_self_induced(
        ::Streamfunction, quad::AbstractQuadrature,
        f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real,
        fit_circle = false,  # ignored
    )
    ℓ₋ = integrate(f, i - 1, quad) do f, j, ζ
        norm(f(j, ζ, Derivative(1)))
    end
    ℓ₊ = integrate(f, i, quad) do f, j, ζ
        norm(f(j, ζ, Derivative(1)))
    end
    # Use local (non-averaged) tangent at point of interest.
    # We want to make sure that the result is exactly tangent to the curve at the node.
    t̂ = f[i, UnitTangent()]
    # Note: the +1 coefficient seems to provide energy conservation to very high accuracy.
    # Tested in particular with leapfrogging rings case.
    β = 2 * prefactor * (log(2 * sqrt(ℓ₋ * ℓ₊) / a) + 1 - Δ)  # note: prefactor = Γ/4π (hence the 2 in front)
    β * t̂
end
