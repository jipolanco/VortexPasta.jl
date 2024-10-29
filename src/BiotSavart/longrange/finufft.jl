# NOTE: this file only defines the types and constants associated to the FINUFFT backends
# (CPU and CUDA), as well as some common definitions for CPU and CUDA backends.
# The actual implementation details are in ext/VortexPastaFINUFFTExt.jl and
# ext/VortexPastaCuFINUFFTExt.jl.
export FINUFFTBackend, CuFINUFFTBackend

using StructArrays: StructArrays, StructArray, StructVector

const FINUFFT_DEFAULT_TOLERANCE = 1e-6
const FINUFFT_DEFAULT_UPSAMPFAC = 1.25  # 1.25 or 2.0

# This includes CPU and CUDA FINUFFT implementations.
abstract type AbstractFINUFFTBackend <: LongRangeBackend end

"""
    FINUFFTBackend <: LongRangeBackend

Compute long-range interactions using the
[FINUFFT.jl](https://github.com/ludvigak/FINUFFT.jl) package.

To use this backend, one first needs to load FINUFFT by doing

    using FINUFFT

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

- `fftw = FFTW.MEASURE` flags passed to the FFTW planner.

Other options described in the [FINUFFT
docs](https://finufft.readthedocs.io/en/latest/opts.html) and not listed above
are also accepted.

"""
struct FINUFFTBackend{KwArgs <: NamedTuple} <: AbstractFINUFFTBackend
    tol :: Float64
    kws :: KwArgs
    # "Private" constructor
    global function _FINUFFTBackend(tol, kws)
        new{typeof(kws)}(tol, kws)
    end
end

"""
    CuFINUFFTBackend <: LongRangeBackend

GPU version of [`FINUFFTBackend`](@ref).

To use this backend, one first needs to load both FINUFFT and CUDA by doing

    using FINUFFT
    using CUDA

Works with Nvidia GPUs only.

!!! compat

    The minimal required version of the FINUFFT libraries is `2.3.0-rc1`.
    In previous versions, cuFINUFFT ignores the `modeord = 1` option which is needed in our
    implementation.

# Optional arguments

The signature of `CuFINUFFTBackend` is:

    CuFINUFFTBackend(; tol = $FINUFFT_DEFAULT_TOLERANCE, kws...)

where all arguments are passed to cuFINUFFT.

Some relevant options are:

- `tol = $FINUFFT_DEFAULT_TOLERANCE` tolerance in NUFFT computations;

- `upsampfac = $FINUFFT_DEFAULT_UPSAMPFAC` upsampling factor. Should be either 1.25 or 2.0;

- `gpu_device::CUDA.CuDevice`: useful if multiple GPUs are available. By default the
  currently active CUDA device is used, i.e. `gpu_device = CUDA.device()`.

See the [cuFINUFFT docs](https://finufft.readthedocs.io/en/latest/c_gpu.html#options-for-gpu-code)
for details and other possible options.

"""
struct CuFINUFFTBackend{
        KwArgs <: NamedTuple,
        Device,  # = CuDevice
        Stream,  # = CuStream
    } <: AbstractFINUFFTBackend
    tol :: Float64
    device :: Device
    stream :: Stream
    kws :: KwArgs
    # "Private" constructor
    global _CuFINUFFTBackend(tol, device, stream, kws) =
        new{typeof(kws), typeof(device), typeof(stream)}(tol, device, stream, kws)
end

## ========================================================================================== ##
## Functions which should be defined by each FINUFFT backend (besides the constructors)

function _finufft_name end
function _finufft_adapt_fourier_vector_field end
function _finufft_unsafe_reshape_vector_to_matrix end
function _finufft_plan_func end
function _finufft_setpts_func! end
function _finufft_exec_func! end
function _finufft_destroy_func! end
function _finufft_sync end

## ========================================================================================== ##
## General definitions for FINUFFTBackend and CuFINUFFTBackend (=> AbstractFINUFFTBackend)

# This is used for both CPU and GPU implementations.
struct FINUFFTCache{
        T,
        CacheCommon <: LongRangeCacheCommon{T},
        Plan,
        ChargeData <: AbstractVector{Complex{T}},
    } <: LongRangeCache
    common :: CacheCommon
    plan_type1 :: Plan  # plan for type-1 NUFFT (physical non-uniform → Fourier uniform)
    plan_type2 :: Plan  # plan for type-2 NUFFT (Fourier uniform → physical non-uniform)
    charge_data :: ChargeData  # used to store charge data used in transforms
end

has_real_to_complex(::AbstractFINUFFTBackend) = false
expected_period(::AbstractFINUFFTBackend) = 2π
# folding_limits(::AbstractFINUFFTBackend) = (-3π, 3π)  # no longer needed since FINUFFT 2.3.0
non_uniform_type(::Type{T}, ::AbstractFINUFFTBackend) where {T <: AbstractFloat} = Complex{T}

# FINUFFT options which should never be modified!
# (Used also in GPU implementation.)
_finufft_options() = (;
    modeord = 1,  # use same mode ordering as FFTW (i.e. k = 0, 1, …, N/2 - 1, -N/2, …, -1)
)

# This should work for CPU and GPU versions.
function init_fourier_vector_field(backend::AbstractFINUFFTBackend, ::Type{T}, Nks::Dims{M}) where {T <: Real, M}
    # Data needs to be in a contiguous array of dimensions (Nks..., M) [usually M = 3].
    device = KA.get_backend(backend)
    data = KA.zeros(device, Complex{T}, (Nks..., M))
    components = ntuple(i -> view(data, :, :, :, i), Val(M))
    uhat = StructArray{Vec3{Complex{T}}}(components)
    @assert uhat isa StructArray{Vec3{Complex{T}}, 3}
    if backend isa FINUFFTBackend  # CPU
        @assert uhat.:1 isa SubArray  # not true on CUDA (which doesn't seem to use the SubArray type)
        @assert parent(uhat.:1) === data  # we use this elsewhere to get the raw data before transforming
    end
    @assert pointer(uhat.:1) == pointer(data)
    uhat
end

# This should work for CPU and GPU versions.
function init_cache_long_ewald(
        pc::ParamsCommon{T},
        params::ParamsLongRange{T, <:AbstractFINUFFTBackend}, args...,
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

    # Make sure that the (cu)finufft_destroy! function is called when the GC frees the
    # plans, so that the memory allocated by the library (in C++) is actually freed.
    # FINUFFT.jl doesn't do this for us...
    destroy! = _finufft_destroy_func!(backend)
    finalizer(destroy!, plan_type1)
    finalizer(destroy!, plan_type2)

    ka_backend = KA.get_backend(backend)
    charge_data = KA.allocate(ka_backend, Complex{T}, 0)  # allocate empty vector
    FINUFFTCache(cache_common, plan_type1, plan_type2, charge_data)
end

# This function is redefined in the VortexPastaThreadPinningExt.jl extension, if
# ThreadPinning.jl is loaded.
# It allows FINUFFT functions to run on unpinned threads, since FINUFFT seems to run much
# slower when ThreadPinning.pinthreads has been used.
# This is only useful on the CPU version; on the GPU version it does nothing.
@inline finufft_unpin_threads(f::F, ::AbstractFINUFFTBackend) where {F} = f()

function _make_finufft_plan_type1(p::AbstractFINUFFTBackend, n_modes::Vector{Int64}, ::Type{T}) where {T}
    opts = _finufft_options()
    type = 1
    iflag = -1  # this sets the FFT convention (we use a - sign for the forward transform)
    ntrans = 3  # number of transforms to compute simultaneously (3 vector components at once)
    _finufft_plan_func(p)(type, n_modes, iflag, ntrans, p.tol; dtype = T, p.kws..., opts...)
end

function _make_finufft_plan_type2(p::AbstractFINUFFTBackend, n_modes::Vector{Int64}, ::Type{T}) where {T}
    opts = _finufft_options()
    type = 2
    iflag = +1  # this sets the FFT convention (we use a + sign for the backward transform)
    ntrans = 3  # number of transforms to compute simultaneously (3 vector components at once)
    _finufft_plan_func(p)(type, n_modes, iflag, ntrans, p.tol; dtype = T, p.kws..., opts...)
end

function _finufft_copy_charges_to_matrix!(
        A::AbstractMatrix{<:Complex},
        qs_in::StructVector{<:Vec3},
    )
    @assert axes(A, 1) == eachindex(qs_in)
    @assert size(A, 2) == 3
    qs = StructArrays.components(qs_in) :: NTuple{3}
    @assert KA.get_backend(A) == KA.get_backend(qs[1])  # avoid implicit host-device copy
    @inbounds for (j, qj) ∈ pairs(qs)
        Aj = @view A[:, j]
        copyto!(Aj, qj)
    end
    A
end

function _finufft_copy_charges_from_matrix!(
        qs_in::StructVector{<:Vec3},
        A::AbstractMatrix{<:Complex},
    )
    @assert axes(A, 1) == eachindex(qs_in)
    @assert size(A, 2) == 3
    qs = StructArrays.components(qs_in) :: NTuple{3}
    @assert KA.get_backend(A) == KA.get_backend(qs[1])  # avoid implicit host-device copy
    @inbounds for (j, qj) ∈ pairs(qs)
        # Note: qj may be a vector of real values, while A is a matrix of complex values.
        # In our case we expect the complex part to be zero.
        Aj = @view A[:, j]
        copyto!(qj, Aj)
    end
    qs_in
end

# Note on synchronisation in the GPU version: on the Julia side, CUDA code runs on the
# stream assigned by CUDA.jl to the currently running Julia task. In particular, if GPU code
# is running asynchronously (e.g. on a @spawn or @async block), it means that it's not
# running on the default Julia task, and therefore it's not running on the same CUDA stream
# as cuFINUFFT code. In other words, we have two different CUDA streams which need to be
# synchronised. Note that KA.synchronize synchronises the current stream from the Julia side,
# while _finufft_sync does the same from the cuFINUFFT side. And we really want the Julia
# stream to be done before calling cuFINUFFT functions.
function transform_to_fourier!(c::FINUFFTCache)
    (; plan_type1, charge_data,) = c
    (; pointdata_d, uhat_d, to_d,) = c.common
    (; points, charges,) = pointdata_d
    backend_lr = backend(c)
    device = KA.get_backend(c)
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points_data = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    Np = length(charges)
    @assert Np == length(points)
    @assert Np > 0
    uhat_data = _finufft_adapt_fourier_vector_field(backend_lr, uhat_d)
    T = eltype(uhat_data)
    @assert T <: Complex
    name_to = _finufft_name(backend_lr)
    KA.synchronize(device)  # make sure point data has been fully written on the GPU
    finufft_unpin_threads(backend_lr) do  # disable ThreadPinning in the CPU version (see VortexPastaThreadPinningExt)
        @timeit to_d "$name_to setpts" begin
            _finufft_setpts_func!(backend_lr)(plan_type1, points_data...)
            # For now we need synchronisation since generally CuFINUFFT code and Julia code
            # run on separate CUDA streams.
            _finufft_sync(backend_lr)  # similar to KA.synchronize, but applies to the CUDA stream attached to cuFINUFFT (irrelevant for CPU case)
        end
        resize!(charge_data, 3 * Np)
        # Note: FINUFFT (CPU version) requires A to be an Array. This means that we can't use Bumper
        # allocators here, as they return some other type of AbstractArray.
        GC.@preserve charge_data begin  # @preserve is only useful on the CPU
            A = _finufft_unsafe_reshape_vector_to_matrix(charge_data, Np, Val(3))
            _finufft_copy_charges_to_matrix!(A, charges)
            _finufft_sync(backend_lr)
            @timeit to_d "$name_to exec" begin
                _finufft_exec_func!(backend_lr)(plan_type1, A, uhat_data)  # execute NUFFT on all components at once
                _finufft_sync(backend_lr)
            end
        end
    end
    _ensure_hermitian_symmetry!(c.common.wavenumbers_d, uhat_d)
    c
end

function _interpolate_to_physical!(output::StructVector, c::FINUFFTCache)
    (; plan_type2, charge_data,) = c
    (; pointdata_d, uhat_d, to_d,) = c.common
    (; points,) = pointdata_d
    backend_lr = backend(c)
    device = KA.get_backend(c)
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points_data = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    Np = length(output)
    @assert Np == length(points)
    uhat_data = _finufft_adapt_fourier_vector_field(backend_lr, uhat_d)
    T = eltype(uhat_data)
    @assert T <: Complex
    name_to = _finufft_name(backend_lr)
    KA.synchronize(device)
    finufft_unpin_threads(backend_lr) do  # disable ThreadPinning in the CPU version (see VortexPastaThreadPinningExt)
        @timeit to_d "$name_to setpts" begin
            _finufft_setpts_func!(backend_lr)(plan_type2, points_data...)
            _finufft_sync(backend_lr)
        end
        resize!(charge_data, 3 * Np)
        GC.@preserve charge_data begin
            A = _finufft_unsafe_reshape_vector_to_matrix(charge_data, Np, Val(3))
            _finufft_sync(backend_lr)
            @timeit to_d "$name_to exec" begin
                _finufft_exec_func!(backend_lr)(plan_type2, uhat_data, A)  # result is computed onto A
                _finufft_sync(backend_lr)
            end
            _finufft_copy_charges_from_matrix!(output, A)
        end
    end
    nothing
end
