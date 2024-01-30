export
    FiniteDiffMethod

@doc raw"""
    FiniteDiffMethod{M} <: DiscretisationMethod{M}
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
    } <: DiscretisationMethod
    interp :: InterpolationMethod

    function FiniteDiffMethod{M}(interp) where {M}
        new{M, typeof(interp)}(interp)
    end

    function FiniteDiffMethod{M, I}() where {M, I}
        new{M, I}(I())
    end
end

continuity(::Type{<:FiniteDiffMethod{M, I}}) where {M, I} = continuity(I)  # returns continuity of interpolation method
npad(::Type{<:FiniteDiffMethod{M}}) where {M} = M

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

struct FiniteDiffCoefs{
        T,
        Method <: FiniteDiffMethod,
        N,  # number of derivatives included (usually 2)
        M,
        Points <: PaddedVector{M, T},
    } <: DiscretisationCoefs{T, Method, N}
    method :: Method

    # Interpolation coefficients associated to the curve.
    # For finite differences, this is simply the node locations.
    cs :: Points

    # Interpolation coefficients associated to first and second derivatives.
    # For finite differences, this is simply the derivatives at nodes.
    cderivs :: NTuple{N, Points}
end

function init_coefficients(
        method::FiniteDiffMethod, cs::V, cderivs::NTuple{N, V},
    ) where {N, V <: PaddedVector}
    _check_coefficients(method, cs, cderivs)
    FiniteDiffCoefs(method, cs, cderivs)
end

allvectors(x::FiniteDiffCoefs) = (x.cs, x.cderivs...)

function compute_coefficients!(
        coefs::FiniteDiffCoefs, Xs::AbstractVector, ts::PaddedVector;
        Xoffset = zero(eltype(Xs)),
    )
    (; method, cs, cderivs,) = coefs
    M = npad(method)
    length(cs) == length(Xs) == length(ts) ||
        throw(DimensionMismatch("incompatible vector dimensions"))
    @assert M == npad(ts) == npad(cs)
    copyto!(cs, Xs)
    pad_periodic!(cs, Xoffset)
    @inbounds for i ∈ eachindex(cs)
        ℓs_i = ntuple(j -> @inbounds(ts[i - M + j] - ts[i - M - 1 + j]), Val(2M))  # = ℓs[(i - M):(i + M - 1)]
        Xs_i = ntuple(j -> @inbounds(cs[i - M - 1 + j]), Val(2M + 1))  # = cs[(i - M):(i + M)]
        coefs_dot = coefs_first_derivative(method, ℓs_i)
        coefs_ddot = coefs_second_derivative(method, ℓs_i)
        cderivs[1][i] = sum(splat(*), zip(coefs_dot, Xs_i))  # = ∑ⱼ c[j] * x⃗[j]
        cderivs[2][i] = sum(splat(*), zip(coefs_ddot, Xs_i))
    end
    # These paddings are needed for Hermite interpolations and stuff like that.
    # (In principle we just need M = 1 for two-point Hermite interpolations.)
    map(pad_periodic!, cderivs)
    nothing
end
