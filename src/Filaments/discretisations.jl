"""
    DiscretisationMethod

Abstract type defining a filament discretisation method.
"""
abstract type DiscretisationMethod end

# This should be implemented by all discretisation methods
continuity(m::DiscretisationMethod) = continuity(typeof(m))
npad(m::DiscretisationMethod) = npad(typeof(m))

"""
    interpolation_method(m::DiscretisationMethod)

Return the interpolation method associated to a [`DiscretisationMethod`](@ref).

For [`FiniteDiffMethod`](@ref) and [`FourierMethod`](@ref), this is usually a
[`HermiteInterpolation`](@ref) (which requires derivatives at the interpolation nodes).
For [`SplineMethod`](@ref)s, since splines don't rely on a separate interpolation scheme,
this simply returns `m` (i.e. the passed `SplineMethod`).
"""
function interpolation_method end

"""
    required_derivatives(m::DiscretisationMethod) -> Int

Return the number of derivatives on interpolation nodes required by the discretisation
method.

This is generally larger than 0 for methods relying on Hermite interpolations (which require
derivatives on interpolation nodes). This is 0 for other methods such as [`SplineMethod`](@ref).
"""
required_derivatives(m::DiscretisationMethod) = required_derivatives(interpolation_method(m))

# Discretisation coefficients associated to a method.
# Here `N` is the number of derivatives that are included in the coefficients.
abstract type DiscretisationCoefs{Method <: DiscretisationMethod, N} end

nderivatives(::Type{<:DiscretisationCoefs{M, N}}) where {M, N} = N
nderivatives(c::DiscretisationCoefs) = nderivatives(typeof(c))

"""
    init_coefficients(method::DiscretisationMethod, ys::AbstractVector, [nderiv::Val]) -> DiscretisationCoefs
    init_coefficients(method::DiscretisationMethod, cs::PaddedVector, cderivs::NTuple) -> DiscretisationCoefs

Initialise interpolation coefficients.

In the first case, the `ys` vector is interpreted as a vector of values at interpolation
nodes, and it will not be modified by the function. Moreover, the `nderivs` argument
indicates the number of derivatives that one wants to be able to compute. By default,
`nderiv` is chosen as the minimum number of derivatives required by the method (see
[`required_derivatives`](@ref)).

In the second case, `cs` is interpreted as a vector of interpolation coefficients, and is
thus stored in the returned structure. Moreover, to compute one or more derivatives, one can
pass `cderivs` vectors which should have the same type and length as `cs`.

Note that some discretisation methods require one or more derivatives, and will fail if
`nderiv = Val(0)` or if `cderivs` is empty. This is the case of methods relying on Hermite
interpolation, which is the default for [`FiniteDiffMethod`](@ref) and
[`FourierMethod`](@ref). One can call [`required_derivatives`](@ref) to know how many
derivatives are required by the discretisation method.
"""
function init_coefficients end

function init_coefficients(
        method::DiscretisationMethod, ys::AbstractVector,
        nderivs::Val = Val(required_derivatives(method)),
    )
    M = npad(method)
    cs = similar_padded(ys, Val(M))
    cderivs = ntuple(_ -> similar(cs), nderivs)
    init_coefficients(method, cs, cderivs)
end

function similar_padded(ys::AbstractVector, ::Val{M}) where {M}
    @assert !(ys isa PaddedVector)
    N = length(ys)
    PaddedVector{M}(similar(ys, N + 2M))
end

similar_padded(ys::PaddedVector{M}, ::Val{M}) where {M} = similar(ys)  # we assume ys has the right padding

function _check_coefficients(
        method::DiscretisationMethod, cs::V, cderivs::NTuple{N, V},
    ) where {N, V <: PaddedVector}
    M = npad(method)
    Nd = required_derivatives(method)
    N ≥ Nd || throw(DimensionMismatch(lazy"$method requires at least $Nd derivatives; got $N derivatives"))
    npad(V) == M ||
        throw(DimensionMismatch(lazy"expected PaddedVector with npad($method) = $M ghost cells in each direction"))
    Np = length(cs)
    all(cp -> length(cp) == Np, cderivs) ||
        throw(DimensionMismatch(lazy"derivative vectors should have length $Np"))
    nothing
end

"""
    compute_coefficients!(coefs::DiscretisationCoefs, [ys::AbstractVector], ts::PaddedVector)

Compute interpolation coefficients of parametric function ``f(t)``.

Here `ys` contains values at interpolation nodes, and `ts` contains the values of the curve
parameter ``t`` at the nodes.

If `ys` is not passed, then the `coefs.cs` vector is expected to have the values at
interpolation points (which will be overwritten with the interpolation coefficients).

This should be called before evaluating an interpolation using [`evaluate`](@ref).
"""
function compute_coefficients! end

function compute_coefficients!(coefs::DiscretisationCoefs, ts::PaddedVector)
    compute_coefficients!(coefs, coefs.cs, ts)
end

"""
    evaluate(coefs::DiscretisationCoefs, ts::PaddedVector, t::Real, [Derivative(0)]; [ileft = nothing])
    evaluate(coefs::DiscretisationCoefs, ts::PaddedVector, i::Int, ζ::Real, [Derivative(0)])

Evaluate interpolation at a given location ``t`` or ``ζ``.

In the first case, ``t`` corresponds to the *global* curve parametrisation.
One can optionally indicate the segment index `i` (such that `ts[i] ≤ t < ts[i + 1]`) usign
the `ileft` keyword argument.

In the second case, ``ζ ∈ [0, 1]`` corresponds to the *local* curve parametrisation, within
the segment `i`.
"""
function evaluate(coefs::DiscretisationCoefs, args...; kws...)
    m = interpolation_method(coefs.method)
    evaluate(m, coefs, args...; kws...)
end

# Check that two sets of coefficients are equal
function Base.:(==)(x::DiscretisationCoefs, y::DiscretisationCoefs)
    x.method === y.method || return false
    us = allvectors(x)
    vs = allvectors(y)
    length(us) === length(vs) || return false
    all(splat(==), zip(us, vs))
end
