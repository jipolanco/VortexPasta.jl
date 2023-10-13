export FourierMethod

using FFTW: FFTW
using LinearAlgebra: mul!, ldiv!

"""
    FourierMethod <: GlobalDiscretisationMethod
    FourierMethod(interpolation = HermiteInterpolation(2))

Describes filaments using Fourier series.

Derivatives at nodes are estimated using FFTs.
To interpolate in-between nodes, a local interpolation method (typically
[`HermiteInterpolation`](@ref)) is used, which is much faster (but less accurate) than
evaluating Fourier series.

Note that using FFTs require the knot locations ``t_i`` to be equidistant.
The default parametrisation used by this method ensures this.
However, this usually works best when the distance between discretisation points is more or
less constant.
"""
struct FourierMethod{
        InterpolationMethod <: LocalInterpolationMethod,
    } <: GlobalDiscretisationMethod
    interp :: InterpolationMethod
    FourierMethod(interp) = new{typeof(interp)}(interp)
    FourierMethod{I}() where {I} = new{I}(I())
end

FourierMethod() = FourierMethod(HermiteInterpolation(2))

# Constant parametrisation with 2π period.
# This shouldn't be modified!! (Assumed in particular when estimating derivatives.)
default_parametrisation(::FourierMethod) = (Xs, i) -> ((i - 1) / length(Xs)) * 2π

interpolation_method(m::FourierMethod) = m.interp

# This is actually the continuity of the (Hermite) interpolation scheme.
continuity(::Type{<:FourierMethod{I}}) where {I} = continuity(I)

npad(::Type{<:FourierMethod}) = 1  # padding needed for Fourier series + Hermite interpolation

Base.show(io::IO, ::FourierMethod{I}) where {I} = print(io, "FourierMethod{", I, "}")

## ================================================================================ ##

mutable struct FourierCoefs{
        Method <: FourierMethod,
        N,  # number of derivatives included (usually 2)
        M,  # padding, needs to be at least 1 for Hermite interpolations
        Points <: PaddedVector{M},
        T <: AbstractFloat,
        Plan,
    } <: DiscretisationCoefs{Method, N}
    const method  :: Method
    const cs      :: Points              # node locations
    const cderivs :: NTuple{N, Points}   # derivatives on nodes
    const ubuf    :: Vector{T}           # FFT input
    const vbuf    :: Vector{Complex{T}}  # FFT output
    plan          :: Plan
end

fftw_flags() = (; flags = FFTW.PRESERVE_INPUT,)

function init_coefficients(method::FourierMethod, Xs::PaddedVector, Nderiv::Val)
    M = npad(Xs)
    @assert M ≥ npad(method)
    cs = similar(Xs)
    cderivs = ntuple(_ -> similar(Xs), Nderiv)
    T = _typeof_number(eltype(Xs))
    ubuf = Vector{T}(undef, length(Xs))
    vbuf = Vector{Complex{T}}(undef, 0)
    plan = FFTW.plan_rfft(ubuf; fftw_flags()...)
    FourierCoefs(method, cs, cderivs, ubuf, vbuf, plan)
end

_typeof_number(::Type{T}) where {T <: Number} = T
_typeof_number(::Type{T}) where {T <: AbstractArray} = _typeof_number(eltype(T))

allvectors(x::FourierCoefs) = (x.cs, x.cderivs...)

## ================================================================================ ##

# Note: `only_derivatives` is not used, it's just there for compatibility with splines.
function _update_coefficients_only!(
        ::FourierMethod, f::ClosedFilament;
        only_derivatives = false,
    )
    (; ts, Xs, Xoffset, coefs,) = f
    (; cs, cderivs, ubuf, vbuf,) = coefs
    M = npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts

    # This function is reused from splines.
    periodise_coordinates!(cs, Xs, ts, Xoffset)

    # TODO optimisations
    # - use StructArrays for coefficients?

    # Compute derivatives in Fourier space.
    N = length(cs)
    ks = FFTW.rfftfreq(N, N)  # assumes parametrisation period T = 2π!!
    Nk = length(ks)
    resize!(ubuf, N)
    resize!(vbuf, Nk)
    if length(coefs.plan) != N
        coefs.plan = FFTW.plan_rfft(ubuf; fftw_flags()...)
    end
    plan = coefs.plan
    @inbounds for i ∈ 1:3
        for j ∈ eachindex(ubuf)
            ubuf[j] = cs[j][i]  # copy input
        end
        mul!(vbuf, plan, ubuf)  # apply FFT
        for n ∈ eachindex(cderivs)
            @. vbuf = im * ks * vbuf
            ldiv!(ubuf, plan, vbuf)  # inverse FFT
            for j ∈ eachindex(ubuf)
                w = cderivs[n]
                w[j] = Base.setindex(w[j], ubuf[j], i)
            end
        end
    end

    # Undo periodisation for the first derivative (higher-order derivatives are not
    # affected).
    cₜ = cderivs[1]
    knot_period = ts[end + 1] - ts[begin]
    @inbounds for j ∈ eachindex(cₜ)
        # This function is reused from splines.
        cₜ[j] = _deperiodise_curve(cₜ[j], Xoffset, knot_period, ts[j], Val(1))
    end

    # These paddings are needed for Hermite interpolations and stuff like that.
    # (In principle we just need M = 1 for two-point Hermite interpolations.)
    map(pad_periodic!, cderivs)

    # Finally, copy coordinates to coefficient vector.
    # This also copies padding (uses copyto! implementation in PaddedArrays)
    copy!(cs, Xs)

    f
end

function _derivative_at_node(
        ::Derivative{n}, ::FourierMethod, f::ClosedFilament, node::AtNode,
    ) where {n}
    (; cs, cderivs,) = f.coefs
    coefs = (cs, cderivs...)
    coefs[n + 1][node.i]
end

# Calls Hermite interpolation functions
_interpolate(m::FourierMethod, args...; kws...) =
    _interpolate(interpolation_method(m), args...; kws...)
