module VortexPastaFINUFFTExt

# Implementation of FINUFFTBackend and CuFINUFFTBackend for computing long-range
# interactions on multiple CPUs (threaded) or on a single CUDA device.
# NOTE: this backend requires at least the cuFINUFFT library v2.3.0(-rc1), which is not
# currently (29/08/2024) included by default with the FINUFFT.jl wrappers. This means that
# one needs to manually compile the FINUFFT libraries and link them to FINUFFT.jl as
# indicated in their README.
# TODO: check version of cuFINUFFT libraries? (how?)

using CUDA
using FINUFFT: FINUFFT
using AbstractFFTs: fftfreq
using FFTW: FFTW
using StructArrays: StructArray
using KernelAbstractions: KernelAbstractions as KA
using StructArrays: StructArrays, StructArray, StructVector
using TimerOutputs: @timeit
using VortexPasta.BiotSavart:
    BiotSavart as BS, Vec3, AbstractFINUFFTBackend, FINUFFTBackend, CuFINUFFTBackend,
    ParamsCommon, ParamsLongRange, LongRangeCacheCommon, LongRangeCache,
    backend

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

BS.has_real_to_complex(::AbstractFINUFFTBackend) = false

BS.expected_period(::AbstractFINUFFTBackend) = 2π
# BS.folding_limits(::AbstractFINUFFTBackend) = (-3π, 3π)  # no longer needed since FINUFFT 2.3.0
BS.non_uniform_type(::Type{T}, ::AbstractFINUFFTBackend) where {T <: AbstractFloat} = Complex{T}

# FINUFFT options which should never be modified!
# (Used also in GPU implementation.)
_finufft_options() = (;
    modeord = 1,  # use same mode ordering as FFTW (i.e. k = 0, 1, …, N/2 - 1, -N/2, …, -1)
)

# This should work for CPU and GPU versions.
function BS.init_fourier_vector_field(backend::AbstractFINUFFTBackend, ::Type{T}, Nks::Dims{M}) where {T <: Real, M}
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
function BS.init_cache_long_ewald(
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
# syncronised. Note that KA.synchronize synchronises the current stream from the Julia side,
# while _finufft_sync does the same from the cuFINUFFT side. And we really want the Julia
# stream to be done before calling cuFINUFFT functions.
function BS.transform_to_fourier!(c::FINUFFTCache)
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
    uhat_data = adapt_fourier_vector_field(backend_lr, uhat_d)
    T = eltype(uhat_data)
    @assert T <: Complex
    name_to = finufft_name(backend_lr)
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
            A = unsafe_reshape_vector_to_matrix(charge_data, Np, Val(3))
            _finufft_copy_charges_to_matrix!(A, charges)
            _finufft_sync(backend_lr)
            @timeit to_d "$name_to exec" begin
                _finufft_exec_func!(backend_lr)(plan_type1, A, uhat_data)  # execute NUFFT on all components at once
                _finufft_sync(backend_lr)
            end
        end
    end
    BS._ensure_hermitian_symmetry!(c.common.wavenumbers_d, uhat_d)
    c
end

function BS._interpolate_to_physical!(output::StructVector, c::FINUFFTCache)
    (; plan_type2, charge_data,) = c
    (; pointdata_d, uhat_d, to_d,) = c.common
    (; points,) = pointdata_d
    backend_lr = backend(c)
    device = KA.get_backend(c)
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    points_data = StructArrays.components(points) :: NTuple{3, <:AbstractVector}
    Np = length(output)
    @assert Np == length(points)
    uhat_data = adapt_fourier_vector_field(backend_lr, uhat_d)
    T = eltype(uhat_data)
    @assert T <: Complex
    name_to = finufft_name(backend_lr)
    KA.synchronize(device)
    finufft_unpin_threads(backend_lr) do  # disable ThreadPinning in the CPU version (see VortexPastaThreadPinningExt)
        @timeit to_d "$name_to setpts" begin
            _finufft_setpts_func!(backend_lr)(plan_type2, points_data...)
            _finufft_sync(backend_lr)
        end
        resize!(charge_data, 3 * Np)
        GC.@preserve charge_data begin
            A = unsafe_reshape_vector_to_matrix(charge_data, Np, Val(3))
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

## ========================================================================================== ##
## Definition of FINUFFTBackend (CPU)

# Define "public" constructor (consistent with the documentation in the BiotSavart module).
function BS.FINUFFTBackend(;
        tol::AbstractFloat = BS.FINUFFT_DEFAULT_TOLERANCE,
        nthreads::Int = Threads.nthreads(),
        upsampfac::Real = BS.FINUFFT_DEFAULT_UPSAMPFAC,
        fftw = FFTW.MEASURE,
        other...,
    )
    kws = (;
        nthreads,
        fftw,
        upsampfac = Float64(upsampfac),
        other...,
    )
    BS._FINUFFTBackend(tol, kws)
end

# This is used for TimerOutputs labels
finufft_name(::FINUFFTBackend) = "FINUFFT"

function Base.show(io::IO, backend::FINUFFTBackend)
    (; tol, kws,) = backend
    (; nthreads, upsampfac,) = kws
    print(io, "FINUFFTBackend(tol = $tol, upsampfac = $upsampfac, nthreads = $nthreads)")
end

# Returns a vector field as a 4D array, as required by FINUFFT.
# Note that this works only on the CPU, and CuFINUFFTBackend uses a different
# implementation.
function adapt_fourier_vector_field(::FINUFFTBackend, uhat::StructArray{Vec3{T}}) where {T}
    ux = uhat.:1  # get first component
    parent(ux) :: Array{T,4}
end

# Equivalent to reshape(v, :, M), but avoids issues when trying to resize `vs` later.
function unsafe_reshape_vector_to_matrix(v::Vector, N, ::Val{M}) where {M}
    @assert length(v) == N * M
    p = pointer(v)
    dims = (N, M)
    unsafe_wrap(Array, p, dims; own = false)
end

_finufft_plan_func(::FINUFFTBackend) = FINUFFT.finufft_makeplan
_finufft_setpts_func!(::FINUFFTBackend) = FINUFFT.finufft_setpts!
_finufft_exec_func!(::FINUFFTBackend) = FINUFFT.finufft_exec!
_finufft_destroy_func!(::FINUFFTBackend) = FINUFFT.finufft_destroy!
_finufft_sync(::FINUFFTBackend) = nothing

## ========================================================================================== ##
## Definition of CuFINUFFTBackend (CUDA)

# Define "public" constructor (consistent with the documentation in the BiotSavart module).
function BS.CuFINUFFTBackend(;
        tol::AbstractFloat = BS.FINUFFT_DEFAULT_TOLERANCE,
        upsampfac = BS.FINUFFT_DEFAULT_UPSAMPFAC,
        gpu_device::CuDevice = CUDA.device(),  # use currently active CUDA device
        # By default use a dedicated CUDA stream to run cuFINUFFT code.
        # Use _finufft_sync to synchronise this stream.
        gpu_stream::CuStream = CUDA.CuStream(),
        other...,
    )
    gpu_stream_ptr = Base.unsafe_convert(Ptr{CUDA.CUstream_st}, gpu_stream)
    kws = (;
        upsampfac = Float64(upsampfac),
        gpu_device_id = CUDA.deviceid(gpu_device),  # extract actual device id (0, 1, ...)
        gpu_stream = gpu_stream_ptr,
        gpu_method = 1,  # method 1 ("GM") seems to be much faster than method 2 ("SM")! (at least with Float64 data)
        gpu_sort = 1,
        gpu_kerevalmeth = 1,
        other...,
    )
    BS._CuFINUFFTBackend(tol, gpu_device, gpu_stream, kws)  # call "private" constructor
end

finufft_name(::CuFINUFFTBackend) = "cuFINUFFT"

# This tells KA to create CUDA kernels, and generic (CPU/GPU) code to create CUDA arrays.
KA.get_backend(::CuFINUFFTBackend) = CUDABackend()

# This longer `show` variant is used when directly printing the backend e.g. in the REPL or
# using @show.
function Base.show(io::IO, ::MIME"text/plain", backend::CuFINUFFTBackend)
    (; tol, device, stream, kws,) = backend
    (; upsampfac,) = kws
    print(io, "CuFINUFFTBackend(tol = $tol, upsampfac = $upsampfac) with:")
    mime = MIME"text/plain"()
    print(io, "\n - CUDA device: ")
    show(io, mime, device)
    print(io, "\n - CUDA stream: ")
    show(io, mime, stream)
end

# This shorter (single-line) variant is used when `print`ing the backend, e.g. as part of a
# larger structure.
function Base.show(io::IO, backend::CuFINUFFTBackend)
    (; tol, kws,) = backend
    (; upsampfac,) = kws
    print(io, "CuFINUFFTBackend(tol = $tol, upsampfac = $upsampfac)")
end

function adapt_fourier_vector_field(::CuFINUFFTBackend, uhat::StructArray{Vec3{T}}) where {T}
    # We can't get the underlying 4D array from the StructArray (unlike in the CPU case), so
    # we cheat and use pointers.
    ux = (uhat.:1) :: CuArray{T,3}  # get first component
    p = pointer(ux)
    dims = size(ux)
    unsafe_wrap(CuArray, p, (dims..., 3)) :: CuArray{T,4}
end

# This works correctly on a variable-size CuVector, unlike the case of CPU Vectors.
# So there's nothing unsafe here!
function unsafe_reshape_vector_to_matrix(v::CuVector, N, ::Val{M}) where {M}
    reshape(v, (N, M))
end

_finufft_plan_func(::CuFINUFFTBackend) = FINUFFT.cufinufft_makeplan
_finufft_setpts_func!(::CuFINUFFTBackend) = FINUFFT.cufinufft_setpts!
_finufft_exec_func!(::CuFINUFFTBackend) = FINUFFT.cufinufft_exec!
_finufft_destroy_func!(::CuFINUFFTBackend) = FINUFFT.cufinufft_destroy!
_finufft_sync(backend::CuFINUFFTBackend) = CUDA.synchronize(backend.stream)  # synchronise CUDA stream running cuFINUFFT code

end
