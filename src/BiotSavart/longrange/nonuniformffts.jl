export NonuniformFFTsBackend, HalfSupport

using NonuniformFFTs: NonuniformFFTs, HalfSupport
using AbstractFFTs: fftfreq

"""
    NonuniformFFTsBackend <: LongRangeBackend

Compute long-range interactions using the
[NonuniformFFTs.jl](https://github.com/jipolanco/NonuniformFFTs.jl) package.

This backend may be faster than other NUFFT-based backends since it allows real valued
non-uniform data, meaning that we can use real-to-complex FFTs to accelerate computations.

Computations are parallelised by default using threads.

# Optional arguments

The signature of `NonuniformFFTsBackend` is:

    NonuniformFFTsBackend(; σ = 1.5, m = HalfSupport(4), kws...)

where all arguments are passed to NonuniformFFTs.jl.

Some relevant options are:

- `σ = 1.5`: upsampling factor, which must be larger than 1. Usual values are between 1.25
  (smaller FFTs, less accurate) and 2.0 (larger FFTs, more accurate). Other values such as 1.5
  (default) also work;

- `m = HalfSupport(4)`: the half-width of the NUFFT kernels. Larger means higher accuracy;

- `fftw_flags = FFTW.MEASURE`: flags passed to the FFTW planner.

The default parameters (`σ = 1.5`, `m = HalfSupport(4)`) correspond to a relative NUFFT
tolerance of ``∼10^{-6}``.
"""
struct NonuniformFFTsBackend{
        HS <: HalfSupport, KwArgs <: NamedTuple,
    } <: LongRangeBackend
    m :: HS
    kws :: KwArgs
    function NonuniformFFTsBackend(;
            σ = 1.5,
            m = HalfSupport(4),
            fftw_flags = FFTW.MEASURE,
            other...,
        )
        kws = (; σ, fftw_flags, other...,)
        hs = to_halfsupport(m)
        new{typeof(hs), typeof(kws)}(hs, kws)
    end
end

to_halfsupport(M::Int) = HalfSupport(M)
to_halfsupport(m::HalfSupport) = m

function Base.show(io::IO, backend::NonuniformFFTsBackend)
    (; m, kws,) = backend
    (; σ,) = kws
    print(io, "NonuniformFFTsBackend(; m = $m, σ = $σ)")
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
    (; m, kws,) = backend
    d = length(Ns)  # dimensionality (usually 3)
    plan = NonuniformFFTs.PlanNUFFT(T, Ns; ntransforms = Val(d), m, kws...)  # plan for real-to-complex transform
    wavenumbers = ntuple(Val(d)) do i
        i == 1 ? rfftfreq(Ns[i], 2π * Ns[i] / Ls[i]) : fftfreq(Ns[i], 2π * Ns[i] / Ls[i])
    end
    cache_common = LongRangeCacheCommon(pc, params, wavenumbers, args...)
    NonuniformFFTsCache(cache_common, plan)
end

function transform_to_fourier!(c::NonuniformFFTsCache)
    (; plan,) = c
    (; pointdata, uhat,) = c.common
    (; points, charges,) = pointdata
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    charges_data = StructArrays.components(charges)
    uhat_data = StructArrays.components(uhat)
    NonuniformFFTs.set_points!(plan, points)
    NonuniformFFTs.exec_type1!(uhat_data, plan, charges_data)  # execute NUFFT on all components at once
    _ensure_hermitian_symmetry!(c.common.wavenumbers, uhat)
    c
end

function interpolate_to_physical!(c::NonuniformFFTsCache)
    (; plan,) = c
    (; pointdata, uhat,) = c.common
    (; points, charges,) = pointdata
    # Interpret StructArrays as tuples of arrays (which is their actual layout).
    charges = StructArrays.components(charges)
    uhat = StructArrays.components(uhat)
    NonuniformFFTs.set_points!(plan, points)
    NonuniformFFTs.exec_type2!(charges, plan, uhat)
    c
end
