# NOTE: this file only defines the types and constants associated to the FINUFFT backends
# (CPU and CUDA). The actual implementation details are in ext/VortexPastaFINUFFTExt.jl.
export FINUFFTBackend, CuFINUFFTBackend

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

To use this backend, one first needs to load FINUFFT by doing

    using FINUFFT

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

function finufft_unpin_threads end
