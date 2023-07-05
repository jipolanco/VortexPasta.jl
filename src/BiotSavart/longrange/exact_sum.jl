export ExactSumBackend

using AbstractFFTs: fftfreq, rfftfreq
using LinearAlgebra: ⋅
using Polyester: @batch  # @threads replacement

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

function init_cache_long_ewald(
        pc::ParamsCommon{T},
        params::ParamsLongRange{<:ExactSumBackend}, timer::TimerOutput,
    ) where {T}
    (; Ls,) = pc
    (; Ns,) = params
    wavenumbers = ntuple(Val(3)) do i
        N, L = Ns[i], Ls[i]
        f = i == 1 ? rfftfreq : fftfreq
        f(N, 2π * N / L)
    end
    cache_common = LongRangeCacheCommon(pc, params, wavenumbers, timer)
    ExactSumCache(cache_common)
end

function reset_fields!(c::ExactSumCache)
    (; uhat,) = c
    fill!(uhat, zero(eltype(uhat)))
    c
end

function set_num_points!(c::ExactSumCache, Np)
    resize!(c.charges, Np)
    c
end

# Add contribution of point charge to Fourier space field `uhat`.
function add_pointcharge!(c::ExactSumCache, X::Vec3, Q::Vec3, i::Int)
    (; uhat, wavenumbers,) = c
    @assert size(uhat) == map(length, wavenumbers)
    inds = CartesianIndices(uhat)
    @inbounds @batch for I ∈ inds
        k⃗ = Vec3(map(getindex, wavenumbers, Tuple(I)))
        uhat[I] += Q * cis(-k⃗ ⋅ X)
    end
    c
end

# Compute Fourier interpolation at X
function add_point!(c::ExactSumCache, X::Vec3, i::Int)
    (; uhat, wavenumbers, charges,) = c
    kxs = first(wavenumbers)
    kx_lims = first(kxs), last(kxs)
    @assert kxs[2] > 0  # only positive half is included (Hermitian symmetry)
    q⃗ = zero(real(eltype(uhat)))
    @inbounds for I ∈ CartesianIndices(uhat)
        k⃗ = Vec3(map(getindex, wavenumbers, Tuple(I)))
        v = uhat[I]
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
    charges[i] = q⃗
    c
end

# Set to zero "asymmetric" modes, to ease comparison with other implementations.
function _ensure_hermitian_symmetry!(c::ExactSumCache, us::AbstractArray)
    _ensure_hermitian_symmetry!(c, Val(ndims(us)), us)
end

# Ensure Hermitian symmetry one dimension at a time.
@inline function _ensure_hermitian_symmetry!(c::ExactSumCache, ::Val{d}, us) where {d}
    N = size(us, d)
    kd = c.wavenumbers[d]
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

# Almost nothing to do, most of the work was already done in `add_pointcharge!`.
# We just zero out some "asymmetric" modes to ease the comparison with other implementations.
function transform_to_fourier!(c::ExactSumCache)
    _ensure_hermitian_symmetry!(c, c.uhat)
    c
end

# Nothing to do, all the work was already done in `add_point!`.
interpolate_to_physical!(c::ExactSumCache) = c
