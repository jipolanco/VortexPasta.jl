import Base: +, *, /, ^, <

"""
    Zero <: Real

Singleton type representing zero.
"""
struct Zero <: Real end

(::Zero) * (::Number) = Zero()
(::Number) * (::Zero) = Zero()
(::Zero) * (::Zero) = Zero()
(::Zero) + (x::Number) = x
(x::Number) + (::Zero) = x
(::Zero) + (::Zero) = Zero()
(::Zero) / (::Number) = Zero()
(::Zero)^(::Integer) = Zero()

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
(::Infinity) * (::Number) = Infinity()
(::Number) * (::Infinity) = Infinity()
(::Infinity) * (::Infinity) = Infinity()
(::Infinity) / (::Number) = Infinity()

(::Number) < (::Infinity) = true
(::Infinity) < (::Real) = false   # needs to be Real to disambiguate
(::Infinity) < (::Infinity) = false

(::Number) / (::Infinity) = Zero()
(::Number) / (::Zero) = Infinity()
