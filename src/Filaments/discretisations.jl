"""
    DiscretisationMethod

Abstract type defining a filament discretisation method.
"""
abstract type DiscretisationMethod end

# This should be implemented by all discretisation methods
continuity(m::DiscretisationMethod) = continuity(typeof(m))

npad(m::DiscretisationMethod) = npad(typeof(m))

"""
    LocalDiscretisationMethod{M} <: DiscretisationMethod

Defines a *local* filament discretisation method.

Local discretisations only require neighbouring node positions to estimate
filament derivatives at a given point.
Finite difference discretisations are a good example of local discretisation.

The parameter `M` is the number of neighbouring nodes that are required on each
side to estimate the derivative.
For instance, `M = 2` for 5-point finite differences (see [`FiniteDiffMethod`](@ref)).
"""
abstract type LocalDiscretisationMethod{M} <: DiscretisationMethod end
npad(::Type{<:LocalDiscretisationMethod{M}}) where {M} = M

"""
    GlobalDiscretisationMethod <: DiscretisationMethod

Defines a *global* filament discretisation method.

Global discretisations require knowing all node positions to estimate filament
derivatives at a given point.

Some examples of this are spline- or Fourier-based filament descriptions.
"""
abstract type GlobalDiscretisationMethod <: DiscretisationMethod end

# Discretisation coefficients associated to a method.
# Here `N` is the number of derivatives that are included in the coefficients.
abstract type DiscretisationCoefs{Method <: DiscretisationMethod, N} end

nderivatives(::Type{<:DiscretisationCoefs{M, N}}) where {M, N} = N
nderivatives(c::DiscretisationCoefs) = nderivatives(typeof(c))

# Check that two sets of coefficients are equal
function Base.:(==)(x::DiscretisationCoefs, y::DiscretisationCoefs)
    x.method === y.method || return false
    us = allvectors(x)
    vs = allvectors(y)
    length(us) === length(vs) || return false
    all(splat(==), zip(us, vs))
end
