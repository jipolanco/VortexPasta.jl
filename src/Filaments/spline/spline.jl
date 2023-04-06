export
    CubicSplineMethod

"""
    CubicSplineMethod <: GlobalDiscretisationMethod

Represents curves using cubic splines.

In the case of closed curves, periodic cubic splines are used.
"""
struct CubicSplineMethod <: GlobalDiscretisationMethod end

init(::Type{ClosedFilament{T}}, N::Integer, ::CubicSplineMethod) where {T} =
    ClosedSplineFilament(N, T)
