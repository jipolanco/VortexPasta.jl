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
"""
abstract type ShortRangeBackend end

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
