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

Base.convert(::Type{T}, ::Zero) where {T <: Real} = zero(T)
Base.convert(::Type{Zero}, ::Zero) = Zero()

Base.:(*)(::Zero, ::Real) = Zero()
Base.:(*)(::Real, ::Zero) = Zero()
Base.:(*)(::Zero, ::Zero) = Zero()
Base.:(+)(x::Real, ::Zero) = x
Base.:(-)(x::Real, ::Zero) = x
Base.:(+)(::Zero, ::Zero) = Zero()
Base.:(/)(::Zero, ::Real) = Zero()
Base.:(÷)(::Zero, ::Real) = Zero()
Base.:(-)(::Zero) = Zero()

Base.exp(::Zero) = true  # in the sense of the multiplicative identity (= 1)

"""
    One <: RealConst <: Real

Singleton type representing the number one.
"""
struct One <: RealConst end

Base.:(*)(::One, x::Real) = x
Base.:(*)(::One, ::Zero) = Zero()

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
Base.:(*)(::Real, ::Infinity) = Infinity()
Base.:(*)(::Infinity, ::Real) = Infinity()
Base.:(+)(::Real, ::Infinity) = Infinity()
Base.:(+)(::Infinity, ::Real) = Infinity()
Base.:(-)(::Infinity, ::Real) = Infinity()
Base.:(-)(::Infinity, ::Zero) = Infinity()
Base.:(+)(::Infinity, ::Infinity) = Infinity()
Base.:(+)(::Infinity, ::Zero) = Infinity()
Base.:(*)(::Infinity, ::Infinity) = Infinity()
Base.:(*)(::One, ::Infinity) = Infinity()
Base.:(/)(::Infinity, ::Real) = Infinity()
Base.:(÷)(::Infinity, ::Real) = Infinity()
Base.:(/)(::Real, ::Infinity) = Zero()

Base.:(<)(::Real, ::Infinity) = true
Base.:(≤)(::Real, ::Infinity) = true
Base.:(>)(::Real, ::Infinity) = false
Base.:(≥)(::Real, ::Infinity) = false
Base.:(<)(::Infinity, ::Real) = false   # needs to be Real to disambiguate
Base.:(<)(::Infinity, ::Infinity) = false
Base.:(>)(::Infinity, ::Infinity) = false

Base.:(/)(::Real, ::Zero) = Infinity()
Base.:(/)(::Infinity, ::Zero) = Infinity()
Base.:(/)(::Zero, ::Infinity) = Zero()

Base.:(/)(::Zero, ::Zero) = error("got Zero() / Zero()")
Base.:(/)(::Infinity, ::Infinity) = error("got Infinity() / Infinity()")
Base.:(*)(::Infinity, ::Zero) = error("got Infinity() * Zero()")
Base.:(*)(::Zero, ::Infinity) = error("got Infinity() * Zero()")

Base.:(≤)(::Vec{N, <:AbstractFloat}, ::Infinity) where {N} = one(Vec{N, Bool})  # true
Base.:(*)(::Zero, ::Vec) = Zero()
Base.:(*)(::One, x::Vec) = x

end
