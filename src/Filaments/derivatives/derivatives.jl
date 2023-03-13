# TODO
# - Implement estimation from local spline interpolation (e.g. 5th order / 5 points).
#   Some tests show that, with "free" boundary conditions, this gives a much
#   better estimation than finite differences with the same number of points.
#   (Note that we can't use natural BCs with odd-order splines...)
#   Not sure about performance though, since spline interpolations need to
#   solve linear systems (small and banded).

"""
    DerivativeEstimationMethod{M}

Abstract type defining a method of estimating curve derivatives at filament nodes.

The parameter `M` is the number of neighbouring nodes that are required on each
side to estimate the derivative.
For instance, `M = 2` for 5-point finite differences (see [`FiniteDiff`](@ref)).
"""
abstract type DerivativeEstimationMethod{M} end

npad(::Type{<:DerivativeEstimationMethod{M}}) where {M} = M
npad(m::DerivativeEstimationMethod) = npad(typeof(m))

include("finitediff.jl")
