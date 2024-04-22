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

A `BACKEND <: ShortRangeBackend` must implement the function:

    init_cache_short(c::ParamsCommon, p::ParamsShortRange{<:BACKEND}, fs::AbstractVector{<:AbstractFilament}, to::TimerOutput)

which should return a [`ShortRangeCache`](@ref).

It may also implement the function [`max_cutoff_distance`](@ref).
"""
abstract type ShortRangeBackend end

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

max_cutoff_distance(backend::ShortRangeBackend, Ls::NTuple) =
    max_cutoff_distance(backend, min(Ls...))

"""
    ShortRangeCache

Abstract type describing the storage of data required to compute short-range interactions.

The [`init_cache_short`](@ref) function returns a concrete instance of a `ShortRangeCache`.

# Interface

## Fields

The following fields must be included in a cache:

- `params :: ParamsShortRange` parameters for short-range computations;

- `to :: TimerOutput` for measuring time spent on different functions.

"""
abstract type ShortRangeCache end
