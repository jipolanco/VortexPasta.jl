module VortexPastaFINUFFTExt

# Implementation of FINUFFTBackend and CuFINUFFTBackend for computing long-range
# interactions on multiple CPUs (threaded) or on a single CUDA device.
# NOTE: this backend requires at least the cuFINUFFT library v2.3.0(-rc1), which is not
# currently (29/08/2024) included by default with the FINUFFT.jl wrappers. This means that
# one needs to manually compile the FINUFFT libraries and link them to FINUFFT.jl as
# indicated in their README.
# TODO: check version of cuFINUFFT libraries? (how?)

using FINUFFT: FINUFFT
using FFTW: FFTW
using StructArrays: StructArray
using VortexPasta.BiotSavart: BiotSavart as BS, Vec3, FINUFFTBackend

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

function Base.show(io::IO, backend::FINUFFTBackend)
    (; tol, kws,) = backend
    (; nthreads, upsampfac,) = kws
    print(io, "FINUFFTBackend(tol = $tol, upsampfac = $upsampfac, nthreads = $nthreads)")
end

# This is used for TimerOutputs labels
BS._finufft_name(::FINUFFTBackend) = "FINUFFT"

# Returns a vector field as a 4D array, as required by FINUFFT.
# Note that this works only on the CPU, and CuFINUFFTBackend uses a different
# implementation.
function BS._finufft_adapt_fourier_vector_field(::FINUFFTBackend, uhat::StructArray{Vec3{T}}) where {T}
    ux = uhat.:1  # get first component
    parent(ux) :: Array{T,4}
end

# Equivalent to reshape(v, :, M), but avoids issues when trying to resize `vs` later.
function BS._finufft_unsafe_reshape_vector_to_matrix(v::Vector, N, ::Val{M}) where {M}
    @assert length(v) == N * M
    p = pointer(v)
    dims = (N, M)
    unsafe_wrap(Array, p, dims; own = false)
end

BS._finufft_plan_func(::FINUFFTBackend) = FINUFFT.finufft_makeplan
BS._finufft_setpts_func!(::FINUFFTBackend) = FINUFFT.finufft_setpts!
BS._finufft_exec_func!(::FINUFFTBackend) = FINUFFT.finufft_exec!
BS._finufft_destroy_func!(::FINUFFTBackend) = FINUFFT.finufft_destroy!
BS._finufft_sync(::FINUFFTBackend) = nothing

end
