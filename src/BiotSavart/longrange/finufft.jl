export FINUFFTBackend

using FINUFFT: FINUFFT, finufft_makeplan, finufft_setpts!, finufft_exec!
using FFTW: FFTW
using AbstractFFTs: fftfreq
using StructArrays: StructArrays, StructVector, StructArray
using StaticArrays: SVector
using LinearAlgebra: ×

const FINUFFT_DEFAULT_TOLERANCE = 1e-8
const FINUFFT_DEFAULT_UPSAMPFAC = 1.25  # 1.25 or 2.0

"""
    FINUFFTBackend <: LongRangeBackend

Compute long-range interactions using the
[FINUFFT.jl](https://github.com/ludvigak/FINUFFT.jl) package.

This package provides a Julia interface to the
[FINUFFT](https://github.com/flatironinstitute/finufft) C++ library,
which enables efficient and accurate computation of non-uniform fast Fourier
transforms (NUFFTs) based on an "exponential of a semi-circle" kernel.

Computations can be parallelised using threads.

# Optional arguments

The signature of `FINUFFTBackend` is:

    FINUFFTBackend(; tol = $FINUFFT_DEFAULT_TOLERANCE, kws...)

where all arguments are passed to FINUFFT.

Some relevant options are:

- `tol = $FINUFFT_DEFAULT_TOLERANCE` tolerance in NUFFT computations;

- `upsampfac = $FINUFFT_DEFAULT_UPSAMPFAC` upsampling factor. Must be either
  1.25 (usually faster) or 2.0 (required to exceed 9 digits of accuracy);

- `nthreads = Threads.nthreads()` number of threads to use. By default, all
  threads available to Julia are used;

- `fftw = FFTW.MEASURE` flags passed to the FFTW planner;

- `chkbnds = false` if `true`, check that non-uniform points are in ``[-3π, 3π)``;

Other options described in the [FINUFFT
docs](https://finufft.readthedocs.io/en/latest/opts.html) and not listed above
are also accepted.

"""
struct FINUFFTBackend{KwArgs <: NamedTuple} <: LongRangeBackend
    tol :: Float64
    kws :: KwArgs
    function FINUFFTBackend(;
            tol::Float64 = FINUFFT_DEFAULT_TOLERANCE,
            nthreads::Int = Threads.nthreads(),
            upsampfac = FINUFFT_DEFAULT_UPSAMPFAC,
            fftw = FFTW.MEASURE,
            chkbnds = false,
            other...,
        )
        kws = (;
            nthreads, fftw, upsampfac = Float64(upsampfac),
            chkbnds = Int(chkbnds), other...,
        )
        new{typeof(kws)}(tol, kws)
    end
end

expected_period(::FINUFFTBackend) = 2π
folding_limits(::FINUFFTBackend) = (-3π, 3π)  # we could even reduce this...

struct FINUFFTCache{
        T,
        CacheCommon <: LongRangeCacheCommon{T},
        Points <: StructVector{Vec3{T}},
        Plan,
    } <: LongRangeCache
    common :: CacheCommon
    plan_type1 :: Plan  # plan for type-1 NUFFT (physical non-uniform → Fourier uniform)
    plan_type2 :: Plan  # plan for type-2 NUFFT (Fourier uniform → physical non-uniform)
    points  :: Points   # non-uniform locations in physical space (3 × [Np])
end

# FINUFFT options which should never be modified!
_finufft_options() = (;
    modeord = 1,  # use same mode ordering as FFTW (i.e. k = 0, 1, …, N/2 - 1, -N/2, …, -1)
)

function _make_finufft_plan_type1(p::FINUFFTBackend, n_modes::Vector{Int64}, ::Type{T}) where {T}
    opts = _finufft_options()
    type = 1
    iflag = -1  # this sets the FFT convention (we use a - sign for the forward transform)
    ntrans = 1  # number of transforms to compute simultaneously (1 vector component at a time)
    finufft_makeplan(type, n_modes, iflag, ntrans, p.tol; dtype = T, p.kws..., opts...)
end

function _make_finufft_plan_type2(p::FINUFFTBackend, n_modes::Vector{Int64}, ::Type{T}) where {T}
    opts = _finufft_options()
    type = 2
    iflag = +1  # this sets the FFT convention (we use a + sign for the backward transform)
    ntrans = 1  # number of transforms to compute simultaneously (1 vector component at a time)
    finufft_makeplan(type, n_modes, iflag, ntrans, p.tol; dtype = T, p.kws..., opts...)
end

function init_cache_long_ewald(
        pc::ParamsCommon{T},
        params::ParamsLongRange{<:FINUFFTBackend}, timer::TimerOutput,
    ) where {T}
    (; Ls,) = pc
    (; backend, Ns,) = params
    n_modes = collect(Int64, Ns)  # type expected by finufft_makeplan
    wavenumbers = map((N, L) -> fftfreq(N, 2π * N / L), Ns, Ls)
    cache_common = LongRangeCacheCommon(pc, params, wavenumbers, timer)
    Nks = map(length, wavenumbers)  # in this case (complex-to-complex transform) this is the same as Ns
    @assert Ns == Nks
    plan_type1 = _make_finufft_plan_type1(backend, n_modes, T)
    plan_type2 = _make_finufft_plan_type2(backend, n_modes, T)
    points = StructVector{Vec3{T}}(undef, 0)
    FINUFFTCache(cache_common, plan_type1, plan_type2, points)
end

function set_num_points!(c::FINUFFTCache, Np)
    resize!(c.points, Np)
    resize!(c.charges, Np)
    c
end

# This is used for type-1 NUFFTs (physical to Fourier).
function add_pointcharge!(c::FINUFFTCache, X::Vec3, Q::Vec3, i::Int)
    @inbounds c.points[i] = X
    @inbounds c.charges[i] = Q
    c
end

# This is used for type-2 NUFFTs (X is an interpolation point).
function add_point!(c::FINUFFTCache, X::Vec3, i::Int)
    @inbounds c.points[i] = X
    c
end

# Set to zero asymmetric modes from complex-to-complex FFT.
function _ensure_hermitian_symmetry!(c::FINUFFTCache, us::Array{<:Complex})
    _ensure_hermitian_symmetry!(c, Val(ndims(us)), us)
end

# Ensure Hermitian symmetry one dimension at a time.
@inline function _ensure_hermitian_symmetry!(c::FINUFFTCache, ::Val{d}, us) where {d}
    N = size(us, d)
    kd = c.wavenumbers[d]
    Δk = kd[2]
    if iseven(N)
        imin = (N ÷ 2) + 1  # asymmetric mode
        @assert -kd[imin] ≈ kd[imin - 1] + Δk
        inds = ntuple(j -> j == d ? imin : Colon(), Val(ndims(us)))
        @inbounds @views us[inds...] .= 0
    end
    _ensure_hermitian_symmetry!(c, Val(d - 1), us)
end

_ensure_hermitian_symmetry!(c::FINUFFTCache, ::Val{0}, us) = us  # we're done, do nothing

function transform_to_fourier!(c::FINUFFTCache)
    (; plan_type1, points, charges, uhat,) = c
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    charges = StructArrays.components(charges)
    uhat = StructArrays.components(uhat)
    finufft_setpts!(plan_type1, points...)
    for (qs, us) ∈ zip(charges, uhat)
        finufft_exec!(plan_type1, qs, us)
        _ensure_hermitian_symmetry!(c, us)
    end
    c
end

function interpolate_to_physical!(c::FINUFFTCache)
    (; plan_type2, points, charges, uhat,) = c
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    charges = StructArrays.components(charges)
    uhat = StructArrays.components(uhat)
    finufft_setpts!(plan_type2, points...)
    for (qs, us) ∈ zip(charges, uhat)
        finufft_exec!(plan_type2, us, qs)
    end
    c
end
