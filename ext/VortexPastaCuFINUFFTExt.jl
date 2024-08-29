module VortexPastaCuFINUFFTExt

# Implementation of CuFINUFFTBackend for computing long-range interactions on a single CUDA
# device.
# NOTE: this backend requires at least the cuFINUFFT library v2.3.0(-rc1), which is not
# currently (29/08/2024) included by default with the FINUFFT.jl wrappers. This means that
# one needs to manually compile the FINUFFT libraries and link them to FINUFFT.jl as
# indicated in their README.
# TODO: check version of cuFINUFFT libraries? (how?)

using CUDA
using FINUFFT
using KernelAbstractions: KernelAbstractions as KA
using StructArrays: StructArrays, StructArray, StructVector
using VortexPasta.BiotSavart: BiotSavart as BS, Vec3, CuFINUFFTBackend

# cuFINUFFT is much slower with oversampling factors different from 2.
const CUFINUFFT_DEFAULT_UPSAMPFAC = 2.0

# This tells KA to create CUDA kernels, and generic (CPU/GPU) code to create CUDA arrays.
KA.get_backend(::CuFINUFFTBackend) = CUDABackend()

# Define "public" constructor (consistent with the documentation in the BiotSavart module).
function BS.CuFINUFFTBackend(;
        tol::AbstractFloat = BS.FINUFFT_DEFAULT_TOLERANCE,
        upsampfac = CUFINUFFT_DEFAULT_UPSAMPFAC,
        gpu_device::CUDA.CuDevice = CUDA.device(),  # use currently active CUDA device
        other...,
    )
    kws = (;
        upsampfac = Float64(upsampfac),
        gpu_device_id = CUDA.deviceid(gpu_device),  # extract actual device id (0, 1, ...)
        other...,
    )
    BS._CuFINUFFTBackend(tol, kws)  # call "private" constructor
end

function Base.show(io::IO, backend::CuFINUFFTBackend)
    (; tol, kws,) = backend
    (; upsampfac,) = kws
    print(io, "CuFINUFFTBackend(tol = $tol, upsampfac = $upsampfac)")
end

function BS.adapt_fourier_vector_field(::CuFINUFFTBackend, uhat::StructArray{Vec3{T}}) where {T}
    # We can't get the underlying 4D array from the StructArray (unlike in the CPU case), so
    # we cheat and use pointers.
    ux = (uhat.:1) :: CuArray{T,3}  # get first component
    p = pointer(ux)
    dims = size(ux)
    unsafe_wrap(CuArray, p, (dims..., 3)) :: CuArray{T,4}
end

BS._finufft_plan_func(::CuFINUFFTBackend) = cufinufft_makeplan
BS._finufft_setpts_func!(::CuFINUFFTBackend) = cufinufft_setpts!
BS._finufft_exec_func!(::CuFINUFFTBackend) = cufinufft_exec!

# This works correctly on a variable-size CuVector, unlike the case of CPU Vectors.
# So there's nothing unsafe here!
function BS.unsafe_reshape_vector_to_matrix(v::CuVector, N, ::Val{M}) where {M}
    reshape(v, (N, M))
end

end
