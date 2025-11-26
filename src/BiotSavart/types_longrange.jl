"""
    LongRangeBackend

Abstract type denoting the backend to use for computing long-range interactions.

# Implemented backends

- [`NonuniformFFTsBackend`](@ref): estimates long-range interactions via the non-uniform fast
  Fourier transform (NUFFT) using the
  [NonuniformFFTs.jl](https://github.com/jipolanco/NonuniformFFTs.jl) package;

- [`ExactSumBackend`](@ref): computes long-range interactions using exact Fourier sums. This
  is really inefficient and should only be used for testing.

# Extended help

## Implementation details

The following functions must be implemented by a `BACKEND <: LongRangeBackend`:

- `init_cache_long_ewald(c::ParamsCommon, p::ParamsLongRange{<:BACKEND}, to::TimerOutput) -> LongRangeCache`,

- [`has_real_to_complex`](@ref),

- [`expected_period`](@ref) (optional),

- [`KernelAbstractions.get_backend`](@ref) (required for GPU-based backends),

- [`KernelAbstractions.device`](@ref) (required for GPU-based backends).

"""
abstract type LongRangeBackend <: AbstractBackend end

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

Currently, this function always returns `true`. It used to return `false` for the
`FINUFFTBackend`, which has been removed.
"""
function has_real_to_complex end

# This is a dummy backend associated to a NullLongRangeCache (meaning that long-range
# computations are disabled).
struct NullLongRangeBackend <: LongRangeBackend end

"""
    expected_period(::LongRangeBackend) -> Union{Nothing, Real}

Domain period expected by the backend.

This is used for rescaling input coordinates to the requirements of the backend.
For instance, NonuniformFFTs.jl assumes a period ``2π``, and therefore coordinates are
rescaled if the input data has a period different from ``2π``.
"""
expected_period(::LongRangeBackend) = nothing

"""
    LongRangeCache

Abstract type describing the storage of data required to compute long-range interactions.

The [`init_cache_long`](@ref) function returns a concrete instance of a `LongRangeCache`
(or `NullLongRangeCache()`, if long-range computations were disabled by setting `α = Zero()`).

# Useful fields

Most useful fields of a `cache::LongRangeCache` are in the `cache.common` field.
In particular, `cache.common` contains the fields:

- `wavenumbers::NTuple{3, AbstractVector}`: Fourier wavenumbers in each direction;

- `uhat::StructArray{Vec3{Complex{T}}, 3}`: a full vector field in Fourier space;

- `pointdata::PointData`: data associated to vector charges applied on non-uniform points.
  These are available in `pointdata.charges` and `pointdata.points`;

Note that, for GPU-based backends, practically all arrays (`uhat` in particular) are GPU arrays,
which don't support all operations of CPU arrays (such as `for` loops).

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

# This is for convenience: doing c.α is equivalent to c.common.α (we do the same for ParamsBiotSavart).
@inline function Base.getproperty(c::LongRangeCache, name::Symbol)
    common = getfield(c, :common)
    if hasproperty(common, name)
        getfield(common, name)
    else
        getfield(c, name)
    end
end

function Base.propertynames(c::LongRangeCache, private::Bool = false)
    (fieldnames(typeof(c))..., propertynames(c.common, private)...)
end

get_parameters(c::LongRangeCache) = c.common.params_all::ParamsBiotSavart
backend(c::LongRangeCache) = backend(get_parameters(c).longrange)
ewald_smoothing_scale(c::LongRangeCache) = get_parameters(c).σ
has_real_to_complex(c::LongRangeCache) = has_real_to_complex(c.common)
KA.get_backend(c::LongRangeCache) = KA.get_backend(backend(c))
KA.device(c::LongRangeCache) = KA.device(backend(c))

add_point_charges!(c::LongRangeCache, fs::AbstractVector{<:AbstractFilament}) =
    add_point_charges!(c.common.pointdata, fs, c.common.params_all)
