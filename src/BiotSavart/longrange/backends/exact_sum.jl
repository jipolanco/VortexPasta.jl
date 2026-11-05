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
        T, N,
        CacheCommon <: LongRangeCacheCommon{T, N},
    } <: LongRangeCache{N}
    common :: CacheCommon
end

has_real_to_complex(::ExactSumBackend) = true

function _init_cache_long_ewald(
        backend::ExactSumBackend,
        params_all::AbstractParams{T},
        params::AbstractParamsLongRange{T}, args...,
    ) where {T}
    @assert backend === params.backend
    pc = params_all.common
    @assert params === params_all.longrange
    (; Ls,) = pc
    (; Ns,) = params
    D = length(Ns)
    wavenumbers = ntuple(Val(D)) do i
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
    kd = c.common.wavenumbers[d]
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

# 3D version
function _transform_to_fourier!(c::ExactSumCache{T, 3}, prefactor::Real) where {T}
    (; uhat, wavenumbers, pointdata,) = c.common
    (; charges,) = pointdata
    xp_src = get_source_points_longrange(pointdata)
    @assert size(uhat) == map(length, wavenumbers)
    fill!(uhat, zero(eltype(uhat)))
    inds = CartesianIndices(uhat)
    @inbounds for i ∈ eachindex(xp_src, charges)
        X = xp_src[i]
        Q = prefactor * charges[i]
        @inbounds Threads.@threads for I ∈ inds
            k⃗ = Vec3(map(getindex, wavenumbers, Tuple(I)))
            uhat[I] += Q * cis(-k⃗ ⋅ X)
        end
    end
    # We zero out some "asymmetric" modes to ease the comparison with other implementations.
    _ensure_hermitian_symmetry!(c, c.common.uhat)
    c
end

# Case of a vector field (either 3D Biot-Savart or 2D velocity field in point vortex)
function _interpolate_to_physical!(callback::F, output::StructVector{<:SVector{N}}, c::ExactSumCache) where {F, N}
    (; uhat, wavenumbers, pointdata,) = c.common
    xp_dest = get_destination_points_longrange(pointdata)
    VecN = SVector{N}
    @assert length(xp_dest) == length(output)
    kxs = first(wavenumbers)
    kx_lims = first(kxs), last(kxs)
    @assert kxs[2] > 0  # only positive half is included (Hermitian symmetry)
    @inbounds Threads.@threads for i in eachindex(xp_dest, output)
        X = xp_dest[i]
        q⃗ = zero(real(eltype(uhat)))
        for I in CartesianIndices(uhat)
            k⃗ = VecN(map(getindex, wavenumbers, Tuple(I)))
            v_orig = uhat[I]::VecN
            v = VecN(@inline callback(Tuple(v_orig), Tuple(I)))  # the callback takes tuples and returns tuples
            z = cis(k⃗ ⋅ X)
            δq⃗ = if k⃗[1] in kx_lims
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
