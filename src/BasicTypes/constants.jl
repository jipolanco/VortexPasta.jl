import Base: +, *, /, ^, <

"""
    Zero <: Real

Singleton type representing zero.
"""
struct Zero <: Real end

(::Zero) * (::Number) = Zero()
(::Number) * (::Zero) = Zero()
(x::Number) + (::Zero) = x
(::Zero) + (::Zero) = Zero()
(::Zero) / (::Number) = Zero()

"""
    Infinity <: Real

Singleton type representing infinity.
"""
struct Infinity <: Real end

"""
    ∞

Alias for Infinity().
"""
const ∞ = Infinity()

# Here we assume that Number > 0, which is the case in our application.
(::Number) * (::Infinity) = Infinity()
(::Infinity) * (::Infinity) = Infinity()
(::Infinity) / (::Number) = Infinity()

(::Infinity) < (::Real) = false   # needs to be Real to disambiguate
(::Infinity) < (::Infinity) = false

(::Number) / (::Zero) = Infinity()
