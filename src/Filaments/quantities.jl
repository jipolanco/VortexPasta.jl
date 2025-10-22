using LinearAlgebra: normalize, norm, ⋅, ×

"""
    GeometricQuantity

Abstract type defining a geometric quantity defined for a filament.

Some available geometric quantities include:

- [`UnitTangent`](@ref),

- [`CurvatureVector`](@ref),

- [`CurvatureScalar`](@ref),

- [`CurvatureBinormal`](@ref),

- [`TorsionScalar`](@ref).

Evaluating geometric quantities works in the same way as evaluating derivatives.
"""
abstract type GeometricQuantity end

Base.broadcastable(q::GeometricQuantity) = Ref(q)

# This changes the order of arguments to make dispatch easier.
@inline function (f::AbstractFilament)(node::AtNode, q::GeometricQuantity)
    _evaluate(q, f, node)
end

@inline function (f::AbstractFilament)(i::Int, ζ::Number, q::GeometricQuantity)
    _evaluate(q, f, i, ζ)
end

@inline function (f::AbstractFilament)(t::Number, q::GeometricQuantity; ileft = nothing)
    _evaluate(q, f, t; ileft)
end

@doc raw"""
    UnitTangent <: GeometricQuantity

Represents the unit tangent vector ``\bm{T}`` at a filament location.

In terms of an arbitrary parametrisation ``\bm{s}(t)`` (where ``t`` is in general different
from the arc length ``ξ``), the unit tangent vector is

```math
\bm{T} = \frac{\bm{s}'}{|\bm{s}'|}
```

where derivatives are with respect to the parameter ``t``.
"""
struct UnitTangent <: GeometricQuantity end

@inline function _evaluate(::UnitTangent, f::AbstractFilament, args...; kws...)
    s′ = f(args..., Derivative(1); kws...)
    normalize(s′)
end

@doc raw"""
    CurvatureVector <: GeometricQuantity

Represents the curvature vector associated to a filament.

In terms of an arbitrary parametrisation ``\bm{s}(t)`` (where ``t`` is in general different
from the arc length ``ξ``), the curvature vector is

```math
\bm{ρ}
= \frac{{|\bm{s}'|}^2 \bm{s}'' - (\bm{s}' ⋅ \bm{s}'') \, \bm{s}'}{|\bm{s}'|^4}
= \frac{\bm{s}' × (\bm{s}'' × \bm{s}')}{|\bm{s}'|^4}
```

where derivatives are with respect to the parameter ``t``.

The curvature vector can be written as ``\bm{ρ} = ρ \bm{N}`` where ``\bm{N}`` is the unit normal
vector and ``ρ`` the scalar curvature (the inverse of the curvature radius).
"""
struct CurvatureVector <: GeometricQuantity end

@inline function _evaluate(::CurvatureVector, f::AbstractFilament, args...; kws...)
    s′ = f(args..., Derivative(1); kws...)
    s″ = f(args..., Derivative(2); kws...)
    s′² = sum(abs2, s′)
    (s′ × (s″ × s′)) ./ (s′² * s′²)
end

@doc raw"""
    CurvatureScalar <: GeometricQuantity

Represents the scalar curvature associated to a filament.

This is simply the norm of [`CurvatureVector`](@ref) and of [`CurvatureBinormal`](@ref).
It is explicitly given by

```math
ρ = \frac{|\bm{s}' × \bm{s}''|}{|\bm{s}'|^3}
```

where derivatives are with respect to the arbitrary parametrisation ``\bm{s}(t)``.
"""
struct CurvatureScalar <: GeometricQuantity end

@inline function _evaluate(::CurvatureScalar, f::AbstractFilament, args...; kws...)
    # It's likely cheaper to compute this as the norm of CurvatureBinormal than as the norm
    # of CurvatureVector.
    b⃗ = _evaluate(CurvatureBinormal(), f, args...; kws...)
    norm(b⃗)
end

@doc raw"""
    CurvatureBinormal <: GeometricQuantity

Represents the scaled binormal vector associated to a filament.

The scaled binormal vector is defined as
``\bm{b} = \bm{T} × ρ⃗ = ρ \, (\bm{T} × \bm{N}) = ρ \bm{B}``,
where ``\bm{B}`` is the (unit) binormal vector and ``ρ`` is the scalar curvature.

It is explicitly given by

```math
\bm{b} = \frac{\bm{s}' × \bm{s}''}{|\bm{s}'|^3}
```

where derivatives are with respect to the arbitrary parametrisation ``\bm{s}(t)``.
"""
struct CurvatureBinormal <: GeometricQuantity end

@inline function _evaluate(::CurvatureBinormal, f::AbstractFilament, args...; kws...)
    s′ = f(args..., Derivative(1); kws...)
    s″ = f(args..., Derivative(2); kws...)
    (s′ × s″) / norm(s′)^3
end

@doc raw"""
    TorsionScalar <: GeometricQuantity

Torsion of a filament.

The torsion ``τ`` describes the variation of the binormal vector along the curve.

It can be obtained as

```math
τ = \frac{(\bm{s}' × \bm{s}'') ⋅ \bm{s}'''}{|\bm{s}' × \bm{s}''|^2}
```

where derivatives are with respect to an arbitrary curve parametrisation.

!!! note

    Because it is obtained from third derivatives, estimating the torsion requires a
    high-order filament discretisation scheme such as [`QuinticSplineMethod`](@ref).

"""
struct TorsionScalar <: GeometricQuantity end

@inline function _evaluate(::TorsionScalar, f::AbstractFilament, args...; kws...)
    s′ = f(args..., Derivative(1); kws...)
    s″ = f(args..., Derivative(2); kws...)
    s‴ = f(args..., Derivative(3); kws...)
    b⃗ = s′ × s″
    (b⃗ ⋅ s‴) / sum(abs2, b⃗)
end
