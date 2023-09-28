export
    FiniteDiffMethod

@doc raw"""
    FiniteDiffMethod{M} <: LocalDiscretisationMethod{M}
    FiniteDiffMethod([M = 2], [interpolation = HermiteInterpolation(M)])

Estimation of curve derivatives at filament nodes using finite differences.

For now, only the case `M = 2` (default) is implemented (4th order / 5-point finite differences),
following the method used by [Baggaley2011](@citet) based on a paper by [Gamet1999](@citet).

`FiniteDiffMethod` also requires specifying an interpolation scheme for evaluating
coordinates and derivatives in-between discretisation points.
By default, Hermite interpolations of continuity ``C^M`` are used.
For the case `M = 2`, that means quintic Hermite interpolations, which match the first two
derivatives estimated by finite differences at the discretisation points.

# References

- [Baggaley2011](@cite) Baggaley & Barenghi, Phys. Rev. B **83**, 134509 (2011)
- [Gamet1999](@cite) Gamet *et al.*, Int. J. Numer. Meth. Fluids **29**, 2 (1999)
"""
struct FiniteDiffMethod{
        M,
        InterpolationMethod <: LocalInterpolationMethod,
    } <: LocalDiscretisationMethod{M}
    interp :: InterpolationMethod

    function FiniteDiffMethod{M}(interp) where {M}
        new{M, typeof(interp)}(interp)
    end
end

continuity(::Type{<:FiniteDiffMethod{M, I}}) where {M, I} = continuity(I)  # returns continuity of interpolation method

@inline FiniteDiffMethod(M::Int, interp = HermiteInterpolation(M)) = FiniteDiffMethod{M}(interp)
@inline FiniteDiffMethod() = FiniteDiffMethod(2)

interpolation_method(m::FiniteDiffMethod) = m.interp

Base.show(io::IO, ::FiniteDiffMethod{M}) where {M} = print(io, "FiniteDiffMethod(", M, ")")

@inline function coefs_first_derivative(::FiniteDiffMethod{2}, ℓ::NTuple{4})
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

@inline function coefs_second_derivative(::FiniteDiffMethod{2}, ℓ::NTuple{4})
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
