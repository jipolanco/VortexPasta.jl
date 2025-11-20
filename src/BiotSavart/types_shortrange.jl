"""
    ShortRangeBackend

Abstract type denoting the backend used for computing short-range interactions.

# Implemented backends

- [`CellListsBackend`](@ref): most efficient when the cutoff radius is much smaller than the
  domain size. Can only be used with periodic boundary conditions.

- [`NaiveShortRangeBackend`](@ref): usually less efficient as it needs to compute distances
  between all pairs of locations.

# Extended help

## Implementation details

The following functions must be implemented by a `BACKEND <: ShortRangeBackend`:

- [`init_cache_short`](@ref),

- [`max_cutoff_distance`](@ref) (optional),

- [`KernelAbstractions.get_backend`](@ref) (required for GPU-based backends),

- [`KernelAbstractions.device`](@ref) (required for GPU-based backends).
"""
abstract type ShortRangeBackend <: AbstractBackend end

"""
    max_cutoff_distance(::ShortRangeBackend, L::Real) -> r
    max_cutoff_distance(::ShortRangeBackend, Ls::NTuple{3, Real}) -> r

Return the maximum cut-off distance `r_cut` allowed by the backend for a given domain period
`L`.

This is usually close to `L/2`, but the actual value can depend on implementation details of
the backend.
For example, the [`CellListsBackend`](@ref) requires a slightly smaller distance, in the range
`L/3 ≤ r_max < L/2` depending on the backend parameters.
"""
function max_cutoff_distance end

# This is the default, can be overridden by specific backends.
max_cutoff_distance(::ShortRangeBackend, L::AbstractFloat) = L / 2

# This is for non-periodic (open) domains.
max_cutoff_distance(::ShortRangeBackend, ::Infinity) = Infinity()

max_cutoff_distance(backend::ShortRangeBackend, Ls::NTuple) =
    max_cutoff_distance(backend, min(Ls...))

"""
    ShortRangeCache

Abstract type describing the storage of data required to compute short-range interactions.

The [`init_cache_short`](@ref) function returns a concrete instance of a `ShortRangeCache`.

# Interface

## Fields

The following fields must be included in a cache:

- `params::ParamsShortRange`: parameters for short-range computations;

- `pointdata::PointData`: updated quadrature locations, charge values and output points (filament nodes);

- `to::TimerOutput`: for measuring time spent on different functions.

"""
abstract type ShortRangeCache end

# This is for convenience: doing c.α is equivalent to c.common.α (we do the same for ParamsBiotSavart).
@inline function Base.getproperty(c::ShortRangeCache, name::Symbol)
    common = getfield(c, :common)
    if hasproperty(common, name)
        getfield(common, name)
    else
        getfield(c, name)
    end
end

function Base.propertynames(c::ShortRangeCache, private::Bool = false)
    (fieldnames(typeof(c))..., propertynames(c.common, private)...)
end

"""
    init_cache_short(pc::ParamsCommon, p::ParamsShortRange{<:ShortRangeBackend}, pointdata::PointData) -> ShortRangeCache

Initialise the cache for the short-range backend defined in `p`.

This should be defined for each [`ShortRangeBackend`](@ref). In general, implementations
should _create a copy_ of `pointdata` to avoid possible aliasing issues, since the
`pointdata` in the signature is usually the one owned by the full [`BiotSavartCache`](@ref).
"""
function init_cache_short end

backend(c::ShortRangeCache) = backend(c.params::ParamsShortRange)
KA.get_backend(c::ShortRangeCache) = KA.get_backend(backend(c))
KA.device(c::ShortRangeCache) = KA.device(backend(c))

"""
    process_point_charges!(cache::ShortRangeCache)

Process list of point charges.

This is useful for short-range backends like [`CellListsBackend`](@ref), which needs to
assign a cell to each point charge before finding nearby pairs.

Must be called after [`add_point_charges!`](@ref) and before computing any short-range quantities
(using [`add_pair_interactions!`](@ref)).
"""
process_point_charges!(::ShortRangeCache) = nothing  # can be overridden by the backend
