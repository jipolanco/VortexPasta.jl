# TODO
# - Implement estimation from local spline interpolation (e.g. 5th order / 5 points).
#   Some tests show that, with "free" boundary conditions, this gives a much
#   better estimation than finite differences with the same number of points.
#   (Note that we can't use natural BCs with odd-order splines...)
#   Not sure about performance though, since spline interpolations need to
#   solve linear systems (small and banded).

"""
    DiscretisationMethod

Abstract type defining a filament discretisation method.
"""
abstract type DiscretisationMethod end

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
