export FINUFFTBackend

using FINUFFT: FINUFFT, finufft_makeplan, finufft_setpts!, finufft_exec!
using FFTW: FFTW
using AbstractFFTs: fftfreq
using StructArrays: StructArrays, StructVector, StructArray
using StaticArrays: SVector
using LinearAlgebra: ×

const FINUFFT_DEFAULT_TOLERANCE = 1e-6
const FINUFFT_DEFAULT_UPSAMPFAC = 1.25  # 1.25 or 2.0

# This includes CPU and CUDA FINUFFT implementations.
abstract type AbstractFINUFFTBackend <: LongRangeBackend end

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
struct FINUFFTBackend{KwArgs <: NamedTuple} <: AbstractFINUFFTBackend
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

function Base.show(io::IO, backend::FINUFFTBackend)
    (; tol, kws,) = backend
    (; nthreads, upsampfac,) = kws
    print(io, "FINUFFTBackend(tol = $tol, upsampfac = $upsampfac, nthreads = $nthreads)")
end

expected_period(::AbstractFINUFFTBackend) = 2π
folding_limits(::AbstractFINUFFTBackend) = (-3π, 3π)  # we could even reduce this...
non_uniform_type(::Type{T}, ::AbstractFINUFFTBackend) where {T <: AbstractFloat} = Complex{T}

function init_fourier_vector_field(::FINUFFTBackend, ::Type{T}, Nks::Dims{M}) where {T <: Real, M}
    # Data needs to be in a contiguous array of dimensions (Nks..., M) [usually M = 3].
    data = Array{Complex{T}}(undef, Nks..., M)
    components = ntuple(i -> view(data, :, :, :, i), Val(M))
    uhat = StructArray{Vec3{Complex{T}}}(components)
    @assert uhat isa StructArray{Vec3{Complex{T}}, 3}
    @assert uhat.:1 isa SubArray
    @assert parent(uhat.:1) === data  # we use this elsewhere to get the raw data before transforming
    uhat
end

struct FINUFFTCache{
        T,
        CacheCommon <: LongRangeCacheCommon{T},
        Plan,
    } <: LongRangeCache
    common :: CacheCommon
    plan_type1 :: Plan  # plan for type-1 NUFFT (physical non-uniform → Fourier uniform)
    plan_type2 :: Plan  # plan for type-2 NUFFT (Fourier uniform → physical non-uniform)
    charge_data :: Vector{Complex{T}}  # used to store charge data used in transforms
end

function init_cache_long_ewald(
        pc::ParamsCommon{T},
        params::ParamsLongRange{T, <:FINUFFTBackend}, args...,
    ) where {T <: AbstractFloat}
    (; Ls,) = pc
    (; backend, Ns,) = params
    n_modes = collect(Int64, Ns)  # type expected by finufft_makeplan
    wavenumbers = map((N, L) -> fftfreq(N, 2π * N / L), Ns, Ls)
    cache_common = LongRangeCacheCommon(pc, params, wavenumbers, args...)
    Nks = map(length, wavenumbers)  # in this case (complex-to-complex transform) this is the same as Ns
    @assert Ns == Nks
    plan_type1 = _make_finufft_plan_type1(backend, n_modes, T)
    plan_type2 = _make_finufft_plan_type2(backend, n_modes, T)
    charge_data = Complex{T}[]
    FINUFFTCache(cache_common, plan_type1, plan_type2, charge_data)
end

# FINUFFT options which should never be modified!
_finufft_options() = (;
    modeord = 1,  # use same mode ordering as FFTW (i.e. k = 0, 1, …, N/2 - 1, -N/2, …, -1)
)

function _make_finufft_plan_type1(p::FINUFFTBackend, n_modes::Vector{Int64}, ::Type{T}) where {T}
    opts = _finufft_options()
    type = 1
    iflag = -1  # this sets the FFT convention (we use a - sign for the forward transform)
    ntrans = 3  # number of transforms to compute simultaneously (3 vector components at once)
    finufft_makeplan(type, n_modes, iflag, ntrans, p.tol; dtype = T, p.kws..., opts...)
end

function _make_finufft_plan_type2(p::FINUFFTBackend, n_modes::Vector{Int64}, ::Type{T}) where {T}
    opts = _finufft_options()
    type = 2
    iflag = +1  # this sets the FFT convention (we use a + sign for the backward transform)
    ntrans = 3  # number of transforms to compute simultaneously (3 vector components at once)
    finufft_makeplan(type, n_modes, iflag, ntrans, p.tol; dtype = T, p.kws..., opts...)
end

function finufft_copy_charges_to_matrix!(
        A::AbstractMatrix{<:Complex},
        qs_in::StructVector{<:Vec3},
    )
    @assert axes(A, 1) == eachindex(qs_in)
    @assert size(A, 2) == 3
    qs = StructArrays.components(qs_in) :: NTuple{3}
    @inbounds for (j, qj) ∈ pairs(qs)
        for (i, q) ∈ pairs(qj)
            A[i, j] = q
        end
    end
    A
end

function finufft_copy_charges_from_matrix!(
        qs_in::StructVector{<:Vec3},
        A::AbstractMatrix{<:Complex},
    )
    @assert axes(A, 1) == eachindex(qs_in)
    @assert size(A, 2) == 3
    qs = StructArrays.components(qs_in) :: NTuple{3}
    @inbounds for (j, qj) ∈ pairs(qs)
        for i ∈ eachindex(qj)
            # Note: qj may be a vector of real values, while A is a matrix of complex values.
            # In our case we expect the complex part to be zero.
            qj[i] = A[i, j]
        end
    end
    qs_in
end

# Equivalent to reshape(v, :, M), but avoids issues when trying to resize `vs` later.
function unsafe_reshape_vector_to_matrix(v::Vector, N, ::Val{M}) where {M}
    @assert length(v) == N * M
    p = pointer(v)
    dims = (N, M)
    unsafe_wrap(Array, p, dims; own = false)
end

function transform_to_fourier!(c::FINUFFTCache)
    (; plan_type1, charge_data,) = c
    (; pointdata_d, uhat_d,) = c.common
    (; points, charges,) = pointdata_d
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points_data = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    Np = length(charges)
    @assert Np == length(points)
    uhat_data = parent(uhat_d.:1)
    T = eltype(uhat_data)
    @assert T <: Complex
    finufft_setpts!(plan_type1, points_data...)
    resize!(charge_data, 3 * Np)
    # Note: FINUFFT requires A to be an Array. This means that we can't use Bumper
    # allocators here, as they return some other type of AbstractArray.
    GC.@preserve charge_data begin
        A = unsafe_reshape_vector_to_matrix(charge_data, Np, Val(3))
        finufft_copy_charges_to_matrix!(A, charges)
        finufft_exec!(plan_type1, A, uhat_data)  # execute NUFFT on all components at once
    end
    _ensure_hermitian_symmetry!(c.common.wavenumbers_d, uhat_d)
    c
end

function interpolate_to_physical!(c::FINUFFTCache)
    (; plan_type2, charge_data,) = c
    (; pointdata_d, uhat_d,) = c.common
    (; points, charges,) = pointdata_d
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points_data = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    Np = length(charges)
    @assert Np == length(points)
    uhat_data = parent(uhat_d.:1)
    T = eltype(uhat_data)
    @assert T <: Complex
    finufft_setpts!(plan_type2, points_data...)
    resize!(charge_data, 3 * Np)
    GC.@preserve charge_data begin
        A = unsafe_reshape_vector_to_matrix(charge_data, Np, Val(3))
        finufft_exec!(plan_type2, uhat_data, A)  # result is computed onto A
        finufft_copy_charges_from_matrix!(charges, A)
    end
    c
end
