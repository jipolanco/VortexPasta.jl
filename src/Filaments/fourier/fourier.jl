export FourierMethod

using FFTW: FFTW
using LinearAlgebra: mul!, ldiv!

"""
    FourierMethod <: DiscretisationMethod
    FourierMethod(interpolation = HermiteInterpolation(2))

Represents closed curves using Fourier series.

Derivatives at nodes are estimated using FFTs.
To interpolate in-between nodes, a local interpolation method (typically
[`HermiteInterpolation`](@ref)) is used, which is much faster (but less accurate) than
evaluating Fourier series.

Note that using FFTs require the knot locations ``t_i`` to be equidistant.
The default parametrisation used by this method ensures this.
However, this usually works best when the distance between discretisation points is more or
less constant.

This method should only be used for simple settings and for verification of filament
derivatives in other methods.
"""
struct FourierMethod{
        InterpolationMethod <: LocalInterpolationMethod,
    } <: DiscretisationMethod
    interp :: InterpolationMethod
    FourierMethod(interp) = new{typeof(interp)}(interp)
    FourierMethod{I}() where {I} = new{I}(I())
end

FourierMethod() = FourierMethod(HermiteInterpolation(2))

default_parametrisation(::FourierMethod) = FourierParametrisation()
interpolation_method(m::FourierMethod) = m.interp

# For practical purposes we set the continuity to a relatively large number.
# Note however that the interpolation method `I` imposes a smaller continuity when doing
# interpolations...
continuity(::Type{<:FourierMethod{I}}) where {I} = 4

npad(::Type{<:FourierMethod}) = 2  # we need at least npad = 1 for Hermite interpolations

Base.show(io::IO, ::FourierMethod{I}) where {I} = print(io, "FourierMethod{", I, "}")

## ================================================================================ ##

mutable struct FourierCoefs{
        T,
        Method <: FourierMethod,
        N,  # number of derivatives included (usually 2)
        M,  # padding, needs to be at least 1 for Hermite interpolations
        Points <: PaddedVector{M, T},
        Tbuf <: AbstractFloat,
        Plan,
    } <: DiscretisationCoefs{T, Method, N}
    const method  :: Method
    const cs      :: Points                 # node locations
    const cderivs :: NTuple{N, Points}      # derivatives on nodes
    const ubuf    :: Vector{Tbuf}           # FFT input
    const vbuf    :: Vector{Complex{Tbuf}}  # FFT output
    plan          :: Plan
end

fftw_flags() = (; flags = FFTW.PRESERVE_INPUT,)

function init_coefficients(
        method::FourierMethod, cs::V, cderivs::NTuple{N, V},
    ) where {N, V <: PaddedVector}
    _check_coefficients(method, cs, cderivs)
    T = _typeof_number(eltype(cs))
    ubuf = Vector{T}(undef, length(cs))
    vbuf = Vector{Complex{T}}(undef, 0)
    plan = FFTW.plan_rfft(ubuf; fftw_flags()...)
    FourierCoefs(method, cs, cderivs, ubuf, vbuf, plan)
end

_typeof_number(::Type{T}) where {T <: Number} = T
_typeof_number(::Type{T}) where {T <: AbstractArray} = _typeof_number(eltype(T))

allvectors(x::FourierCoefs) = (x.cs, x.cderivs...)

function compute_coefficients!(
        coefs::FourierCoefs, Xs::AbstractVector, ts::PaddedVector;
        Xoffset = zero(eltype(Xs)),
        only_derivatives = false,  # unused for now
    )
    (; cs, cderivs, ubuf, vbuf,) = coefs
    M = npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts
    T = eltype(Xs)
    is_scalar = T <: Number  # otherwise, T is usually SVector
    components = eachindex(first(Xs))

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
    @inbounds for i ∈ components
        for j ∈ eachindex(ubuf)
            ubuf[j] = cs[j][i]  # copy input
        end
        mul!(vbuf, plan, ubuf)  # apply FFT
        for n ∈ eachindex(cderivs)
            @. vbuf = im * ks * vbuf
            ldiv!(ubuf, plan, vbuf)  # inverse FFT
            for j ∈ eachindex(ubuf)
                w = cderivs[n]
                w[j] = if is_scalar
                    ubuf[j]
                else
                    Base.setindex(w[j], ubuf[j], i)
                end
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
    foreach(pad_periodic!, cderivs)

    # Finally, copy coordinates to coefficient vector.
    copyto!(cs, Xs)
    pad_periodic!(cs, Xoffset)  # apply periodic padding

    coefs
end

## ================================================================================ ##

# Note: `only_derivatives` is not used, it's just there for compatibility with splines.
function _update_coefficients_only!(
        ::FourierMethod, f::ClosedFilament;
        only_derivatives = false,
    )
    (; ts, Xs, Xoffset, coefs,) = f
    parametrisation(f) === FourierParametrisation() || error("FourierMethod requires using FourierParametrisation")
    compute_coefficients!(coefs, Xs, ts; Xoffset)
    f
end

function _derivative_at_node(
        ::Derivative{n}, ::FourierMethod, f::ClosedFilament, node::AtNode,
    ) where {n}
    (; cs, cderivs,) = f.coefs
    n > length(cderivs) && throw(ArgumentError(
        lazy"derivatives of order $n are not enabled. Initialise a filament with nderivs ≥ $n, and make sure the discretisation scheme has the required continuity degree."
    ))
    coefs = (cs, cderivs...)
    coefs[n + 1][node.i]
end

# Calls Hermite interpolation functions
_interpolate(m::FourierMethod, args...; kws...) =
    _interpolate(interpolation_method(m), args...; kws...)

## ================================================================================ ##
## Refinement

# The refinement of Fourier filaments is very simple: we multiply or divide the number of
# nodes by a power of 2, according to whether the criterium tells us to overall add or
# remove nodes.
# In the case where nodes are added, we do this using Fourier "interpolation", by
# increasing the number of coefficients in Fourier space (i.e. padding with zeroes).
function _refine!(method::FourierMethod, f::AbstractFilament, crit) :: NTuple{2,Int}
    @assert method === discretisation_method(f)

    # Determine where to add or remove nodes.
    cache = _nodes_to_refine!(f, crit)
    (; inds_add, inds_rem,) = cache
    n_add = length(inds_add)
    n_rem = length(inds_rem)
    iszero(n_add + n_rem) && return (n_add, n_rem)

    # In FourierMethod, we just care about the overall change in the number of nodes.
    n_mod = n_add - n_rem
    if n_mod == 0
        return (zero(n_add), zero(n_rem))
    elseif n_mod < 0
        return _fourier_remove_nodes!(f, n_mod)
    else
        return _fourier_insert_nodes!(f, n_mod)
    end
end

function _fourier_remove_nodes!(f::AbstractFilament, n_mod)
    @assert n_mod < 0
    N = length(f)
    N_wanted = N + n_mod
    α = ceil(UInt, N // N_wanted) - 1  # ≥ 1
    # This is the definition of Base.top_set_bit, see also ?Base.top_set_bit.
    # The actual number of nodes will be N / 2ⁿ.
    n = 8 * sizeof(α) - leading_zeros(α)
    @assert n ≥ 1
    # @show N_wanted, N, α, n, 2^n
    n_add = 0
    n_rem = 0
    while n > 0
        # Remove one out of two nodes.
        inds_to_remove = eachindex(f)[begin:2:end]
        for i ∈ reverse(inds_to_remove)
            remove_node!(f, i)
            n_rem += 1
        end
        update_after_changing_nodes!(f; removed = true)
        n -= 1
    end
    n_add, n_rem
end

function _fourier_insert_nodes!(f::AbstractFilament, n_mod)
    @assert n_mod > 0
    N = length(f)
    N_wanted = N + n_mod
    α = floor(UInt, N_wanted // N)  # ≥ 1
    # The actual number of nodes will be N * 2ⁿ.
    n = 8 * sizeof(α) - leading_zeros(α)
    @assert n ≥ 1
    β = 1 << n  # = 2ⁿ
    n_add = (β - 1) * N  # = (2ⁿ - 1) * N
    n_rem = 0
    M = β * N
    Mk = (M >> 1) + 1  # number of Fourier modes (case of real-to-complex transform)
    Nk = (N >> 1) + 1
    Xs = nodes(f)
    resize!(Xs, M)
    (; cs, ubuf, vbuf,) = f.coefs
    sizehint!(ubuf, M)
    sizehint!(vbuf, Mk)
    resize!(ubuf, M)
    renormalisation = M / N  # renormalisation of Fourier coefficients
    plan_old = f.coefs.plan
    plan_new = FFTW.plan_rfft(ubuf; fftw_flags()...)
    f.coefs.plan = plan_new
    @inbounds for i ∈ 1:3
        resize!(ubuf, N)
        resize!(vbuf, Nk)
        for j ∈ eachindex(ubuf)
            ubuf[j] = cs[j][i]  # copy input
        end
        mul!(vbuf, plan_old, ubuf)  # apply FFT
        nlast = lastindex(vbuf)
        resize!(ubuf, M)
        resize!(vbuf, Mk)
        for j ∈ firstindex(vbuf):nlast
            vbuf[j] *= renormalisation
        end
        for j ∈ (nlast + 1):lastindex(vbuf)
            vbuf[j] = 0  # set all added entries to zero
        end
        ldiv!(ubuf, plan_new, vbuf)  # new locations in physical space
        for j ∈ eachindex(Xs, ubuf)
            Xs[j] = Base.setindex(Xs[j], ubuf[j], i)  # sets i-th coordinate (i ∈ 1:3)
        end
    end
    update_after_changing_nodes!(f; removed = false)
    n_add, n_rem
end
