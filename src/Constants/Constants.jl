"""
    Constants

Module defining singleton types representing constants (just `Infinity()` right now).

For now, these are mainly used internally in non-periodic simulations.

This module may dissapear in the future.
"""
module Constants

export Infinity

abstract type RealConst end

"""
    Infinity <: RealConst

Singleton type representing infinity.
"""
struct Infinity <: RealConst end

Base.:(/)(::Infinity, ::Real) = Infinity()
Base.:(*)(::Infinity, ::Infinity) = Infinity()

end
