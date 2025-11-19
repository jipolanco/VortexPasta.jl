export NonuniformFFTsBackend, HalfSupport

using NonuniformFFTs: NonuniformFFTs, HalfSupport
using AbstractFFTs: rfftfreq, fftfreq
using FFTW: FFTW

"""
    NonuniformFFTsBackend <: LongRangeBackend

Compute long-range interactions using the
[NonuniformFFTs.jl](https://github.com/jipolanco/NonuniformFFTs.jl) package.

This backend may be faster than other NUFFT-based backends since it allows real valued
non-uniform data, meaning that we can use real-to-complex FFTs to accelerate computations.

Transforms can be performed either on the CPU (parallelised with threads, default) or on a
single GPU (in principle any kind of GPU should work, but only CUDA has been tested).
This must be set via the first positional argument (see below).

# Optional arguments

The signature of `NonuniformFFTsBackend` is:

    NonuniformFFTsBackend([ka_backend = CPU()]; device = 1, σ = 1.5, m = HalfSupport(4), kws...)

where all arguments besides `device` are passed to NonuniformFFTs.jl.

## Using a GPU

Transforms are run on all available CPUs by default.
To use a GPU, pass the corresponding [KernelAbstractions.jl
backend](https://juliagpu.github.io/KernelAbstractions.jl/stable/#Supported-backends) as the
only positional argument.
For example, to use a CUDA device:

    using CUDA
    backend_long = NonuniformFFTsBackend(CUDABackend(); device = 1, kwargs...)

On AMD GPUs the following should work:

    using AMDGPU
    backend_long = NonuniformFFTsBackend(ROCBackend(); device = 1, kwargs...)

If running on a machine with multiple GPU devices, one may use the `device` keyword argument
to choose the device where long-range computations will be performed. This should be a value
in `1:ndevices`. When using KernelAbstractions.jl, the number of available devices can be
obtained using `KA.ndevices(ka_backend)`. By default the first device (`device = 1`) is used.

## Keyword arguments

Some relevant keyword arguments which are passed to NonuniformFFTs.jl are:

- `σ = 1.5`: upsampling factor, which must be larger than 1. Usual values are between 1.25
  (smaller FFTs, less accurate) and 2.0 (larger FFTs, more accurate). Other values such as 1.5
  (default) also work;

- `m = HalfSupport(4)`: the half-width of the NUFFT kernels. Larger means higher accuracy;

- `use_atomics = Threads.nthreads() > 4`: whether to use atomics (instead of locks) in
  spreading on the CPU (ignored on GPU devices). See NonuniformFFTs docs for details.
  By default we enable this for relatively large number of threads, in which case it can be
  beneficial;

- `fftw_flags = FFTW.MEASURE`: flags passed to the FFTW planner (ignored on GPU devices).

The default parameters (`σ = 1.5`, `m = HalfSupport(4)`) correspond to a relative NUFFT
tolerance of ``∼10^{-6}``.

See [the NonuniformFFTs.jl docs](https://jipolanco.github.io/NonuniformFFTs.jl/stable/API/#NonuniformFFTs.PlanNUFFT)
for a full list of possible keyword arguments.

## Effect of parameters on accuracy

The following table roughly relates accuracy (in number of digits) and NUFFT parameters, as
detailed in [Polanco2025](@citet):

| Precision digits |  NUFFT ``m`` |  NUFFT ``σ`` |  Ewald ``β`` |
| :--------------: | :----------: | :----------: | :----------: |
|         3        |      2       |     1.5      |     2.0      |
|         4        |      3       |     1.5      |     2.5      |
|         6        |      4       |     1.5      |     3.5      |
|         8        |      5       |     1.5      |     4.0      |
|        10        |      6       |     1.5      |     4.5      |
|        12        |      7       |     1.5      |     5.0      |
|        14        |      8       |     1.5      |     5.5      |

The last column is the associated value of the accuracy parameter ``β`` in Ewald's method as
formulated in [Polanco2025](@citet). Once one has set ``β`` and Ewald's splitting parameter
``α`` (an inverse lengthscale), the cut-offs in physical and Fourier space are ``r_{\\text{cut}} = β / α``
and ``k_{\\text{max}} = 2βα``. In this formulation, ``β`` controls the method accuracy while
``α`` is tuned to maximise performance.

"""
struct NonuniformFFTsBackend{
        HS <: HalfSupport, OversamplingFactor <: Real,
        BackendKA <: KA.Backend,
        KwArgs <: NamedTuple,
    } <: LongRangeBackend
    m :: HS
    σ :: OversamplingFactor
    ka_backend :: BackendKA
    ka_device  :: Int
    kws :: KwArgs
    function NonuniformFFTsBackend(
            ka_backend::KA.Backend = ka_default_cpu_backend();
            device::Integer = 1,
            σ = 1.5,
            m = HalfSupport(4),
            fftw_flags = FFTW.MEASURE,
            use_atomics = Threads.nthreads() > 4,
            other...,
        )
        # Pass the chosen KA backend to NonuniformFFTs, except if the backend is a PseudoGPU
        # (used in testing only). Actually passing a PseudoGPU to NonuniformFFTs might work,
        # but would need to be tested...
        backend = ka_backend isa PseudoGPU ? ka_default_cpu_backend() : ka_backend
        kws = (; backend, fftw_flags, use_atomics, other...,)
        hs = to_halfsupport(m)
        KA.device!(backend, device)  # this will fail if `device` is an invalid device id
        new{typeof(hs), typeof(σ), typeof(ka_backend), typeof(kws)}(hs, σ, ka_backend, device, kws)
    end
end

KA.get_backend(backend::NonuniformFFTsBackend) = backend.ka_backend
KA.device(backend::NonuniformFFTsBackend) = backend.ka_device

has_real_to_complex(::NonuniformFFTsBackend) = true

to_halfsupport(M::Int) = HalfSupport(M)
to_halfsupport(m::HalfSupport) = m

oversampling_factor(backend::NonuniformFFTsBackend) = backend.σ
half_support(backend::NonuniformFFTsBackend) = half_support(backend.m)  # returns half support as an integer
half_support(::HalfSupport{M}) where {M} = M

function Base.show(io::IO, backend::NonuniformFFTsBackend)
    (; ka_backend, ka_device, m, σ,) = backend
    print(io, "NonuniformFFTsBackend($ka_backend; device = $ka_device, m = $m, σ = $σ)")
end

expected_period(::NonuniformFFTsBackend) = 2π

# This is not needed since folding is done by NonuniformFFTs anyway:
# folding_limits(::NonuniformFFTsBackend) = (0, 2π)

struct NonuniformFFTsCache{
        T,
        Backend <: NonuniformFFTsBackend,
        CacheCommon <: LongRangeCacheCommon{T},
        Plan,
    } <: LongRangeCache
    backend :: Backend
    common :: CacheCommon
    plan :: Plan  # plan for NUFFTs in both directions
end

function init_cache_long_ewald(
        params_all::ParamsBiotSavart{T},
        params::ParamsLongRange{T, <:NonuniformFFTsBackend}, args...,
    ) where {T}
    pc = params_all.common
    @assert params === params_all.longrange
    (; Ls,) = pc
    (; backend, Ns,) = params
    (; m, σ, kws, ka_backend, ka_device,) = backend
    KA.device!(ka_backend, ka_device)  # change the device if needed
    d = length(Ns)  # dimensionality (usually 3)
    plan = NonuniformFFTs.PlanNUFFT(T, Ns; ntransforms = Val(d), m, σ, kws...)  # plan for real-to-complex transform
    wavenumbers = ntuple(Val(d)) do i
        freq = T(2π * Ns[i] / Ls[i])
        i == 1 ? rfftfreq(Ns[i], freq) : fftfreq(Ns[i], freq)
    end
    cache_common = LongRangeCacheCommon(params_all, wavenumbers, args...)
    NonuniformFFTsCache(backend, cache_common, plan)
end

function transform_to_fourier!(c::NonuniformFFTsCache, prefactor::Real)
    (; backend, plan,) = c
    (; pointdata, uhat,) = c.common
    (; points, charges,) = pointdata
    (; ka_backend, ka_device,) = backend
    # Make sure we're already running on the wanted device (e.g. GPU 2).
    # Usually we call this function right after having defined data (points + charges) on
    # this device, and thus we can expect that we have already selected the right device.
    @assert KA.device(ka_backend) == ka_device
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    charges_data = StructArrays.components(charges)
    uhat_data = StructArrays.components(uhat)
    # Note: we could apply prefactor either in non-uniform (physical) or in uniform (Fourier) space.
    # For now we choose to do it in uniform space, not sure if it's better though.
    callbacks = NonuniformFFTs.NUFFTCallbacks(uniform = @inline((ω̂, idx) -> ω̂ .* prefactor))  # use a callback to apply prefactor
    NonuniformFFTs.set_points!(plan, points)
    NonuniformFFTs.exec_type1!(uhat_data, plan, charges_data; callbacks)  # execute NUFFT on all components at once
    _ensure_hermitian_symmetry!(c.common.wavenumbers, uhat)
    c
end

# Note: the callback must have the signature (û::NTuple{3}, idx::NTuple{3,Int}).
function _interpolate_to_physical!(callback_uniform::F, output::StructVector, c::NonuniformFFTsCache) where {F <: Function}
    (; backend, plan,) = c
    (; pointdata, uhat,) = c.common
    (; nodes,) = pointdata
    (; ka_backend, ka_device,) = backend
    # Make sure we're already running on the wanted device (e.g. GPU 2).
    # Usually we call this function right after having defined data (uhat) on this device,
    # and thus we can expect that we have already selected the right device.
    @assert KA.device(ka_backend) == ka_device
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    charges = StructArrays.components(output)
    uhat_data = StructArrays.components(uhat)
    callbacks = NonuniformFFTs.NUFFTCallbacks(uniform = callback_uniform)
    NonuniformFFTs.set_points!(plan, nodes)
    NonuniformFFTs.exec_type2!(charges, plan, uhat_data; callbacks)
    nothing
end
