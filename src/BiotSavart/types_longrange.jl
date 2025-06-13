"""
    LongRangeBackend

Abstract type denoting the backend to use for computing long-range interactions.

# Implemented backends

- [`NonuniformFFTsBackend`](@ref): estimates long-range interactions via the non-uniform fast
  Fourier transform (NUFFT) using the
  [NonuniformFFTs.jl](https://github.com/jipolanco/NonuniformFFTs.jl) package;

- [`FINUFFTBackend`](@ref): estimates long-range interactions via the NUFFT using the
  [FINUFFT](https://github.com/flatironinstitute/finufft) library;

- [`CuFINUFFTBackend`](@ref): estimates long-range interactions via the NUFFT using the CUDA
  implementation of the [FINUFFT](https://github.com/flatironinstitute/finufft) library.
  CUDA.jl must be loaded before using this backend (CUDA devices only);

- [`ExactSumBackend`](@ref): computes long-range interactions using exact Fourier sums. This
  is really inefficient and should only be used for testing.

# Extended help

## Implementation details

The following functions must be implemented by a `BACKEND <: LongRangeBackend`:

- `init_cache_long_ewald(c::ParamsCommon, p::ParamsLongRange{<:BACKEND}, to::TimerOutput) -> LongRangeCache`.

- [`has_real_to_complex`](@ref),

- [`expected_period`](@ref) (optional),

- [`folding_limits`](@ref) (optional),

- [`KernelAbstractions.get_backend`](@ref) (required for GPU-based backends).

"""
abstract type LongRangeBackend end

"""
    has_real_to_complex(::LongRangeBackend) -> Bool
    has_real_to_complex(::ParamsLongRange) -> Bool
    has_real_to_complex(::LongRangeCacheCommon) -> Bool
    has_real_to_complex(::LongRangeCache) -> Bool

Check whether the backend performs real-to-complex (fast) Fourier transforms.

If `true`, it means that input non-uniform data in physical space can be real-valued, and
that uniform data in Fourier space only contains half the total number of modes along the
first dimension to account for Hermitian symmetry.

This function is useful in particular for:

- knowing which kind of non-uniform data (vorticities) one must give to the backend;
- knowing how to interpret Fourier-space data, e.g. to compute Fourier spectra.

This function returns `false` for backends such as [`FINUFFTBackend`](@ref), as these
require complex input data and don't take advantage of Hermitian symmetry.
"""
function has_real_to_complex end

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

# Useful fields

Most useful fields of a `cache::LongRangeCache` are in the `cache.common` field.
In particular, `cache.common` contains the fields:

- `wavenumbers_d::NTuple{3, AbstractVector}`: Fourier wavenumbers in each direction;
- `uhat_d::StructArray{Vec3{Complex{T}}, 3}`: a full vector field in Fourier space;
- `pointdata_d::PointData`: data associated to vector charges applied on non-uniform points.
  These are available in `pointdata_d.charges` and `pointdata_d.points`;
- `ewald_prefactor::Real`: the quantity ``Γ / V`` where ``V`` is the volume of a periodic
  cell.

The `_d` suffixes means that data is on the computing device associated to the long-range
backend (i.e. on the GPU for GPU-based backends).

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

get_parameters(c::LongRangeCache) = c.common.params_all::ParamsBiotSavart
backend(c::LongRangeCache) = backend(get_parameters(c).longrange)
ewald_smoothing_scale(c::LongRangeCache) = get_parameters(c).σ
has_real_to_complex(c::LongRangeCache) = has_real_to_complex(c.common)
KA.get_backend(c::LongRangeCache) = KA.get_backend(backend(c))

# TODO: this is not optimised for GPU backends
function add_point_charges!(c::LongRangeCache, fs::AbstractVector{<:AbstractFilament})
    (; quad,) = c.common.params
    add_point_charges!(c.common.pointdata_d, fs, quad)
end
