using LinearAlgebra: normalize, norm, ⋅, ×

"""
    GeometricQuantity

Abstract type defining a geometric quantity defined for a filament.

Some available geometric quantities include:

- [`UnitTangent`](@ref),

- [`CurvatureVector`](@ref),

- [`CurvatureScalar`](@ref),

- [`CurvatureBinormal`](@ref).

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

"""
    UnitTangent <: GeometricQuantity

Represents the unit tangent vector ``t̂`` at a filament location.
"""
struct UnitTangent <: GeometricQuantity end

@inline function _evaluate(::UnitTangent, f::AbstractFilament, args...; kws...)
    ṡ = f(args..., Derivative(1); kws...)
    normalize(ṡ)
end

"""
    CurvatureVector <: GeometricQuantity

Represents the curvature vector associated to a filament.

The curvature vector is given by ``\\bm{ρ} = ρ n̂`` where ``n̂`` is the unit normal
vector and ``ρ`` the scalar curvature (the inverse of the curvature radius).
"""
struct CurvatureVector <: GeometricQuantity end

@inline function _evaluate(::CurvatureVector, f::AbstractFilament, args...; kws...)
    ṡ = f(args..., Derivative(1); kws...)
    s̈ = f(args..., Derivative(2); kws...)
    ṡ² = sum(abs2, ṡ)
    (ṡ² * s̈ - (s̈ ⋅ ṡ) * ṡ) ./ (ṡ² * ṡ²)
end

"""
    CurvatureScalar <: GeometricQuantity

Represents the scalar curvature associated to a filament.

This is simply the norm of [`CurvatureVector`](@ref).
"""
struct CurvatureScalar <: GeometricQuantity end

@inline function _evaluate(::CurvatureScalar, f::AbstractFilament, args...; kws...)
    ρ⃗ = _evaluate(CurvatureVector(), f, args...; kws...)
    norm(ρ⃗)
end

"""
    CurvatureBinormal <: GeometricQuantity

Represents the scaled binormal vector associated to a filament.

The scaled binormal vector is defined as ``\\bm{b} = t̂ × ρ⃗ = ρ \\, (t̂ × n̂) = ρ \\, b̂``,
where ``b̂`` is the (unit) binormal vector and ``ρ`` is the scalar curvature.
"""
struct CurvatureBinormal <: GeometricQuantity end

@inline function _evaluate(::CurvatureBinormal, f::AbstractFilament, args...; kws...)
    ṡ = f(args..., Derivative(1); kws...)
    s̈ = f(args..., Derivative(2); kws...)
    (ṡ × s̈) / norm(ṡ)^3
end
