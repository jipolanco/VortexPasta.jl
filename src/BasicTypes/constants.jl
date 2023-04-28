import Base: +, *, /, ^, <

struct Zero <: Real end

(::Zero) * (::Number) = Zero()
(::Number) * (::Zero) = Zero()
(::Zero) * (::Zero) = Zero()
(::Zero) + (x::Number) = x
(x::Number) + (::Zero) = x
(::Zero) + (::Zero) = Zero()
(::Zero) / (::Number) = Zero()
(::Zero)^(::Integer) = Zero()

struct Infinity <: Real end

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
