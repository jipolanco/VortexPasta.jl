import Base: +, -, *, /, ^, <, ÷

abstract type RealConst <: Real end

"""
    Zero <: RealConst <: Real

Singleton type representing zero.
"""
struct Zero <: RealConst end

Base.convert(::Type{T}, ::Zero) where {T <: Number} = zero(T)
Base.convert(::Type{Zero}, ::Zero) = Zero()

(::Zero) * (::Number) = Zero()
(::Number) * (::Zero) = Zero()
(x::Number) + (::Zero) = x
(x::Number) - (::Zero) = x
(::Zero) + (::Zero) = Zero()
(::Zero) / (::Number) = Zero()
(::Zero) ÷ (::Number) = Zero()

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

Base.isinf(::Infinity) = true

# Here we assume that Number > 0, which is the case in our application.
(::Number) * (::Infinity) = Infinity()
(::Infinity) * (::Infinity) = Infinity()
(::Infinity) / (::Number) = Infinity()

(::Real) > (::Infinity) = false
(::Infinity) < (::Real) = false   # needs to be Real to disambiguate
(::Infinity) < (::Infinity) = false

(::Number) / (::Zero) = Infinity()
