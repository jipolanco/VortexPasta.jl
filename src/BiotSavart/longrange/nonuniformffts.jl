export NonuniformFFTsBackend, HalfSupport

using NonuniformFFTs: NonuniformFFTs, HalfSupport
using AbstractFFTs: fftfreq

"""
    NonuniformFFTsBackend <: LongRangeBackend

Compute long-range interactions using the
[NonuniformFFTs.jl](https://github.com/jipolanco/NonuniformFFTs.jl) package.

This backend may be faster than other NUFFT-based backends since it allows real valued
non-uniform data, meaning that we can use real-to-complex FFTs to accelerate computations.

Transforms can be performed either on the CPU (parallelised with threads, default) or on a
single GPU (in principle any kind of GPU should work, but only CUDA has been tested).
This must be set via the `device` argiment (see below).

# Optional arguments

The signature of `NonuniformFFTsBackend` is:

    NonuniformFFTsBackend([device = CPU()]; σ = 1.5, m = HalfSupport(4), kws...)

where all arguments are passed to NonuniformFFTs.jl.

## Using a GPU

Transforms are run on all available CPUs by default.
To use a GPU, pass the corresponding [KernelAbstractions.jl
backend](https://juliagpu.github.io/KernelAbstractions.jl/stable/#Supported-backends) as the
only positional argument.
For example, to use a CUDA device:

    using CUDA
    backend_long = NonuniformFFTsBackend(CUDABackend(); kwargs...)

On AMD GPUs, the following should hopefully work (not tested):

    using AMDGPU
    backend_long = NonuniformFFTsBackend(ROCBackend(); kwargs...)

## Keyword arguments

Some relevant keyword arguments are:

- `σ = 1.5`: upsampling factor, which must be larger than 1. Usual values are between 1.25
  (smaller FFTs, less accurate) and 2.0 (larger FFTs, more accurate). Other values such as 1.5
  (default) also work;

- `m = HalfSupport(4)`: the half-width of the NUFFT kernels. Larger means higher accuracy;

- `fftw_flags = FFTW.MEASURE`: flags passed to the FFTW planner (ignored on GPU devices).

The default parameters (`σ = 1.5`, `m = HalfSupport(4)`) correspond to a relative NUFFT
tolerance of ``∼10^{-6}``.

See [the NonuniformFFTs.jl docs](https://jipolanco.github.io/NonuniformFFTs.jl/stable/API/#NonuniformFFTs.PlanNUFFT)
for a full list of possible keyword arguments.
"""
struct NonuniformFFTsBackend{
        HS <: HalfSupport, OversamplingFactor <: Real,
        Device <: KA.Backend,
        KwArgs <: NamedTuple,
    } <: LongRangeBackend
    m :: HS
    σ :: OversamplingFactor
    device :: Device
    kws :: KwArgs
    function NonuniformFFTsBackend(
            device::KA.Backend = ka_default_cpu_backend();
            σ = 1.5,
            m = HalfSupport(4),
            fftw_flags = FFTW.MEASURE,
            other...,
        )
        # Pass the chosen device to NonuniformFFTs, except if the device is a PseudoGPU
        # (used in testing only). Actually passing a PseudoGPU to NonuniformFFTs might work,
        # but would need to be tested...
        backend = device isa PseudoGPU ? ka_default_cpu_backend() : device
        kws = (; backend, fftw_flags, other...,)
        hs = to_halfsupport(m)
        new{typeof(hs), typeof(σ), typeof(device), typeof(kws)}(hs, σ, device, kws)
    end
end

KA.get_backend(backend::NonuniformFFTsBackend) = backend.device

has_real_to_complex(::NonuniformFFTsBackend) = true

to_halfsupport(M::Int) = HalfSupport(M)
to_halfsupport(m::HalfSupport) = m

oversampling_factor(backend::NonuniformFFTsBackend) = backend.σ
half_support(backend::NonuniformFFTsBackend) = half_support(backend.m)  # returns half support as an integer
half_support(::HalfSupport{M}) where {M} = M

function Base.show(io::IO, backend::NonuniformFFTsBackend)
    (; device, m, σ,) = backend
    print(io, "NonuniformFFTsBackend($device; m = $m, σ = $σ)")
end

expected_period(::NonuniformFFTsBackend) = 2π

# This is not needed since folding is done by NonuniformFFTs anyway:
# folding_limits(::NonuniformFFTsBackend) = (0, 2π)

struct NonuniformFFTsCache{
        T,
        CacheCommon <: LongRangeCacheCommon{T},
        Plan,
    } <: LongRangeCache
    common :: CacheCommon
    plan :: Plan  # plan for NUFFTs in both directions
end

function init_cache_long_ewald(
        pc::ParamsCommon{T},
        params::ParamsLongRange{T, <:NonuniformFFTsBackend}, args...,
    ) where {T}
    (; Ls,) = pc
    (; backend, Ns,) = params
    (; m, σ, kws,) = backend
    d = length(Ns)  # dimensionality (usually 3)
    plan = NonuniformFFTs.PlanNUFFT(T, Ns; ntransforms = Val(d), m, σ, kws...)  # plan for real-to-complex transform
    wavenumbers = ntuple(Val(d)) do i
        i == 1 ? rfftfreq(Ns[i], 2π * Ns[i] / Ls[i]) : fftfreq(Ns[i], 2π * Ns[i] / Ls[i])
    end
    cache_common = LongRangeCacheCommon(pc, params, wavenumbers, args...)
    NonuniformFFTsCache(cache_common, plan)
end

function transform_to_fourier!(c::NonuniformFFTsCache)
    (; plan,) = c
    (; pointdata_d, uhat_d,) = c.common
    (; points, charges,) = pointdata_d
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    charges_data = StructArrays.components(charges)
    uhat_data = StructArrays.components(uhat_d)
    NonuniformFFTs.set_points!(plan, points)
    NonuniformFFTs.exec_type1!(uhat_data, plan, charges_data)  # execute NUFFT on all components at once
    _ensure_hermitian_symmetry!(c.common.wavenumbers_d, uhat_d)
    c
end

function _interpolate_to_physical!(output::StructVector, c::NonuniformFFTsCache)
    (; plan,) = c
    (; pointdata_d, uhat_d,) = c.common
    (; points,) = pointdata_d
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    charges = StructArrays.components(output)
    uhat_data = StructArrays.components(uhat_d)
    NonuniformFFTs.set_points!(plan, points)
    NonuniformFFTs.exec_type2!(charges, plan, uhat_data)
    nothing
end
