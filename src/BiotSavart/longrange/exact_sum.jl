export ExactSumBackend

using AbstractFFTs: fftfreq, rfftfreq
using LinearAlgebra: ⋅

"""
    ExactSumBackend <: LongRangeBackend

Compute long-range interactions "exactly" (up to filament discretisation
errors) using sums of Fourier series across non-uniform points.

This should only be used for testing, as it is very slow and scales very badly
with the number of non-uniform points.
"""
struct ExactSumBackend <: LongRangeBackend end

struct ExactSumCache{
        CacheCommon <: LongRangeCacheCommon,
    } <: LongRangeCache
    common :: CacheCommon
end

has_real_to_complex(::ExactSumBackend) = true

function init_cache_long_ewald(
        params_all::ParamsBiotSavart{T},
        params::ParamsLongRange{T, <:ExactSumBackend}, args...,
    ) where {T}
    pc = params_all.common
    @assert params === params_all.longrange
    (; Ls,) = pc
    (; Ns,) = params
    wavenumbers = ntuple(Val(3)) do i
        N, L = Ns[i], Ls[i]
        f = i == 1 ? rfftfreq : fftfreq
        f(N, T(2π * N / L))
    end
    cache_common = LongRangeCacheCommon(params_all, wavenumbers, args...)
    ExactSumCache(cache_common)
end

# Set to zero "asymmetric" modes, to ease comparison with other implementations.
function _ensure_hermitian_symmetry!(c::ExactSumCache, us::AbstractArray)
    _ensure_hermitian_symmetry!(c, Val(ndims(us)), us)
end

# Ensure Hermitian symmetry one dimension at a time.
@inline function _ensure_hermitian_symmetry!(c::ExactSumCache, ::Val{d}, us) where {d}
    N = size(us, d)
    kd = c.common.wavenumbers_d[d]
    Δk = kd[2]
    imin = if kd[end] > 0  # real-to-complex dimension (rfftfreq)
        @assert d == 1
        N
    elseif iseven(N)
        imin_ = (N ÷ 2) + 1  # asymmetric mode
        @assert -kd[imin_] ≈ kd[imin_ - 1] + Δk
        imin_
    else
        0
    end
    if imin > 0
        inds = ntuple(j -> j == d ? imin : Colon(), Val(ndims(us)))
        vs = view(us, inds...)
        fill!(vs, zero(eltype(vs)))
    end
    _ensure_hermitian_symmetry!(c, Val(d - 1), us)
end

_ensure_hermitian_symmetry!(::ExactSumCache, ::Val{0}, us) = us  # we're done, do nothing

function transform_to_fourier!(c::ExactSumCache, prefactor::Real)
    (; uhat_d, wavenumbers_d, pointdata_d,) = c.common
    (; points, charges,) = pointdata_d
    @assert size(uhat_d) == map(length, wavenumbers_d)
    fill!(uhat_d, zero(eltype(uhat_d)))
    inds = CartesianIndices(uhat_d)
    @inbounds for i ∈ eachindex(points, charges)
        X = points[i]
        Q = prefactor * charges[i]
        @inbounds Threads.@threads for I ∈ inds
            k⃗ = Vec3(map(getindex, wavenumbers_d, Tuple(I)))
            uhat_d[I] += Q * cis(-k⃗ ⋅ X)
        end
    end
    # We zero out some "asymmetric" modes to ease the comparison with other implementations.
    _ensure_hermitian_symmetry!(c, c.common.uhat_d)
    c
end

function _interpolate_to_physical!(callback::F, output::StructVector, c::ExactSumCache) where {F}
    (; uhat_d, wavenumbers_d, pointdata_d,) = c.common
    (; points,) = pointdata_d
    @assert length(points) == length(output)
    kxs = first(wavenumbers_d)
    kx_lims = first(kxs), last(kxs)
    @assert kxs[2] > 0  # only positive half is included (Hermitian symmetry)
    @inbounds Threads.@threads for i ∈ eachindex(points, output)
        X = points[i]
        q⃗ = zero(real(eltype(uhat_d)))
        for I ∈ CartesianIndices(uhat_d)
            k⃗ = Vec3(map(getindex, wavenumbers_d, Tuple(I)))
            v_orig = uhat_d[I]::Vec3
            v = @inline callback(v_orig, I)
            z = cis(k⃗ ⋅ X)
            δq⃗ = if k⃗[1] ∈ kx_lims
                # Note: the imaginary part will cancel out with -k⃗ (also computed)
                real(v) * real(z) - imag(v) * imag(z)  # = real(v * z)
            else
                # Apply Hermitian symmetry
                2 * (
                    real(v) * real(z) -
                    imag(v) * imag(z)
                ) # = real(v * z + conj(v) * conj(z))
            end
            q⃗ = q⃗ + δq⃗
        end
        output[i] = q⃗
    end
    nothing
end
