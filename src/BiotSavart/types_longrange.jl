"""
    LongRangeBackend

Abstract type denoting the backend to use for computing long-range interactions.

# Implemented backends

- [`NonuniformFFTsBackend`](@ref): estimates long-range interactions via the non-uniform fast
  Fourier transform (NUFFT) using the
  [NonuniformFFTs.jl](https://github.com/jipolanco/NonuniformFFTs.jl) package;

- [`FINUFFTBackend`](@ref): estimates long-range interactions via the NUFFT using the
  [FINUFFT](https://github.com/flatironinstitute/finufft) library;

- [`ExactSumBackend`](@ref): computes long-range interactions using exact Fourier sums. This
  is really inefficient and should only be used for testing.

# Extended help

## Implementation details

The following functions must be implemented by a `BACKEND <: LongRangeBackend`:

- `init_cache_long_ewald(c::ParamsCommon, p::ParamsLongRange{<:BACKEND}, to::TimerOutput) -> LongRangeCache`.

- [`expected_period`](@ref) (optional),

- [`folding_limits`](@ref) (optional).

"""
abstract type LongRangeBackend end

# This is a dummy backend associated to a NullLongRangeCache (meaning that long-range
# computations are disabled).
struct NullLongRangeBackend <: LongRangeBackend end

"""
    expected_period(::LongRangeBackend) -> Union{Nothing, Real}

Domain period expected by the backend.

This is used for rescaling input coordinates to the requirements of the backend.
For instance, FINUFFT assumes a period ``2π``, and therefore coordinates are
rescaled if the input data has a period different from ``2π``.
"""
expected_period(::LongRangeBackend) = nothing

"""
    folding_limits(::LongRangeBackend) -> Union{Nothing, NTuple{2, Real}}

Domain limits required by the backend.

This is used for folding input coordinates so that they are within the limits
expected by the backend.
For instance, FINUFFT requires coordinates to be in the ``[-3π, 3π]`` interval.

Note that, if a backend defines `folding_limits`, then it must also define
[`expected_period`](@ref).
"""
folding_limits(::LongRangeBackend) = nothing

"""
    LongRangeCache

Abstract type describing the storage of data required to compute long-range interactions.

The [`init_cache_long`](@ref) function returns a concrete instance of a `LongRangeCache`
(or `NullLongRangeCache()`, if long-range computations were disabled by setting `α = Zero()`).

# Extended help

## Implementation details

### Fields

All caches must include a `common <: LongRangeCacheCommon` field which contains common
definitions for all backends.

### Functions

The following functions must be implemented by a cache:

- [`transform_to_fourier!`](@ref),

- [`interpolate_to_physical!`](@ref).

"""
abstract type LongRangeCache end

backend(c::LongRangeCache) = backend(c.common.params)

function add_point_charges!(c::LongRangeCache, fs::AbstractVector{<:AbstractFilament})
    (; quad,) = c.common.params
    add_point_charges!(c.common.pointdata, fs, quad)
end
