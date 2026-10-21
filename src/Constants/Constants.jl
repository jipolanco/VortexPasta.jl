"""
    Constants

Module defining singleton types representing constants (just `Infinity()` right now).

For now, these are mainly used internally in non-periodic simulations.

This module may dissapear in the future.
"""
module Constants

export Infinity

"""
    Infinity <: Real

Singleton type representing (positive) infinity.
"""
struct Infinity <: Real end

Base.isinf(::Infinity) = true
Base.:(/)(::Infinity, ::Union{Integer, AbstractFloat}) = Infinity()  # we could put Real, but this leads to ambiguity with ForwardDiff.Dual
Base.:(*)(::Infinity, ::Infinity) = Infinity()

end
