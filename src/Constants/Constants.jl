"""
    Constants

Module defining singleton types representing constants (0, infinity, …).

For now, these are mainly useful for preparing non-periodic simulations.
"""
module Constants

using SIMD: SIMD, Vec  # for operations mixing constants and SIMD vectors

export Zero, One, Infinity, ∞

abstract type RealConst <: Real end

"""
    Zero <: RealConst <: Real

Singleton type representing zero.
"""
struct Zero <: RealConst end

Base.convert(::Type{T}, ::Zero) where {T <: Number} = zero(T)
Base.convert(::Type{Zero}, ::Zero) = Zero()

Base.:(*)(::Zero, ::Number) = Zero()
Base.:(*)(::Number, ::Zero) = Zero()
Base.:(*)(::Zero, ::Zero) = Zero()
Base.:(+)(x::Number, ::Zero) = x
Base.:(-)(x::Number, ::Zero) = x
Base.:(+)(::Zero, ::Zero) = Zero()
Base.:(/)(::Zero, ::Number) = Zero()
Base.:(÷)(::Zero, ::Number) = Zero()
Base.:(-)(::Zero) = Zero()

Base.exp(::Zero) = true  # in the sense of the multiplicative identity (= 1)

"""
    One <: RealConst <: Real

Singleton type representing the number one.
"""
struct One <: RealConst end

Base.:(*)(::One, x::Number) = x

"""
    Infinity <: RealConst <: Real

Singleton type representing infinity.
"""
struct Infinity <: RealConst end

"""
    ∞

Alias for Infinity().
"""
const ∞ = Infinity()

Base.convert(::Type{T}, ::Infinity) where {T <: AbstractFloat} = T(Inf)

Base.isinf(::Infinity) = true

# Here we assume that Number > 0, which is the case in our application.
Base.:(*)(::Number, ::Infinity) = Infinity()
Base.:(*)(::Infinity, ::Number) = Infinity()
Base.:(+)(::Number, ::Infinity) = Infinity()
Base.:(+)(::Infinity, ::Number) = Infinity()
Base.:(-)(::Infinity, ::Number) = Infinity()
Base.:(+)(::Infinity, ::Infinity) = Infinity()
Base.:(*)(::Infinity, ::Infinity) = Infinity()
Base.:(/)(::Infinity, ::Number) = Infinity()
Base.:(÷)(::Infinity, ::Number) = Infinity()
Base.:(/)(::Number, ::Infinity) = Zero()

Base.:(<)(::Real, ::Infinity) = true
Base.:(≤)(::Real, ::Infinity) = true
Base.:(>)(::Real, ::Infinity) = false
Base.:(≥)(::Real, ::Infinity) = false
Base.:(<)(::Infinity, ::Real) = false   # needs to be Real to disambiguate
Base.:(<)(::Infinity, ::Infinity) = false
Base.:(>)(::Infinity, ::Infinity) = false

Base.:(/)(::Number, ::Zero) = Infinity()

Base.:(≤)(::Vec{N, <:AbstractFloat}, ::Infinity) where {N} = one(Vec{N, Bool})  # true
Base.:(*)(::Zero, ::Vec) = Zero()
Base.:(*)(::One, x::Vec) = x

end
