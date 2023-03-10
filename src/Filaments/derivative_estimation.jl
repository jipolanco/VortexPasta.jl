export
    FiniteDiff

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

@doc raw"""
    FiniteDiff{M} <: DerivativeEstimationMethod{M}

Estimation of curve derivatives at filament nodes using finite differences.

For now, only the case `M = 2` is implemented (4th order / 5-point finite differences),
following the method used by Baggaley & Barenghi[^1] based on a paper by Gamet et al.[^2].

[^1]: A. W. Baggaley & C. F. Barenghi, [Phys. Rev. B **83**, 134509 (2011)](http://dx.doi.org/10.1103/PhysRevB.83.134509).

[^2]: L. Gamet, F. Ducros, F. Nicoud & T. Poinsot, [Int. J. Numer. Meth. Fluids **29**, 2 (1999)](http://dx.doi.org/10.1002/(SICI)1097-0363(19990130)29:2<159::AID-FLD781>3.0.CO;2-9).
"""
struct FiniteDiff{M} <: DerivativeEstimationMethod{M} end

FiniteDiff(M::Int) = FiniteDiff{M}()

Base.show(io::IO, ::FiniteDiff{M}) where {M} = print(io, "FiniteDiff(", M, ")")

@inline function coefs_first_derivative(::FiniteDiff{2}, ℓ::NTuple{4})
    A = ℓ[1] * ℓ[2] * (ℓ[3] + ℓ[4]) / (
        (ℓ[1])
      * (ℓ[1] + ℓ[2])
      * (ℓ[1] + ℓ[2] + ℓ[3])
      * (ℓ[1] + ℓ[2] + ℓ[3] + ℓ[4])
    )
    B = -(
        (ℓ[1] + ℓ[2]) * ℓ[3] * (ℓ[3] + ℓ[4])
    ) / (
        ℓ[1] * ℓ[2] * (ℓ[2] + ℓ[3]) * (ℓ[2] + ℓ[3] + ℓ[4])
    )
    D = (
        (ℓ[1] + ℓ[2]) * ℓ[2] * (ℓ[3] + ℓ[4])
    ) / (
        ℓ[3] * ℓ[4] * (ℓ[2] + ℓ[3]) * (ℓ[1] + ℓ[2] + ℓ[3])
    )
    E = -ℓ[3] * ℓ[2] * (ℓ[2] + ℓ[1]) / (
        (ℓ[4])
      * (ℓ[4] + ℓ[3])
      * (ℓ[4] + ℓ[3] + ℓ[2])
      * (ℓ[4] + ℓ[3] + ℓ[2] + ℓ[1])
    )
    C = -(A + B + D + E)
    (A, B, C, D, E)
end

@inline function coefs_second_derivative(::FiniteDiff{2}, ℓ::NTuple{4})
    A = 2 * (
        ℓ[3] * (-2 * ℓ[2] + ℓ[3])
        +
        ℓ[4] * (-ℓ[2] + ℓ[3])
    ) / (
        (ℓ[1])
      * (ℓ[1] + ℓ[2])
      * (ℓ[1] + ℓ[2] + ℓ[3])
      * (ℓ[1] + ℓ[2] + ℓ[3] + ℓ[4])
    )
    B = 2 * (
        ℓ[3] * (2 * ℓ[1] + 2 * ℓ[2] - ℓ[3])
        +
        ℓ[4] * (ℓ[1] + ℓ[2] - ℓ[3])
    ) / (
        ℓ[1] * ℓ[2] * (ℓ[2] + ℓ[3]) * (ℓ[2] + ℓ[3] + ℓ[4])
    )
    D = 2 * (
        ℓ[1] * (-ℓ[2] + ℓ[3] + ℓ[4])
        + ℓ[2] * (-ℓ[2] + 2 * ℓ[3] + 2 * ℓ[4])
    ) / (
        ℓ[3] * ℓ[4] * (ℓ[2] + ℓ[3]) * (ℓ[1] + ℓ[2] + ℓ[3])
    )
    E = 2 * (
        ℓ[2] * (ℓ[1] + ℓ[2])
        -
        ℓ[3] * (ℓ[1] + 2 * ℓ[2])
    ) / (
        (ℓ[4])
      * (ℓ[4] + ℓ[3])
      * (ℓ[4] + ℓ[3] + ℓ[2])
      * (ℓ[4] + ℓ[3] + ℓ[2] + ℓ[1])
    )
    C = -(A + B + D + E)
    (A, B, C, D, E)
end
