using ..FindNearbySegments:
    NearbySegmentFinder,
    NaiveSegmentFinder,
    CellListSegmentFinder

abstract type AbstractReconnectionCache end

distance(c::AbstractReconnectionCache) = distance(criterion(c))

struct ReconnectionCache{
        Criterion <: ReconnectionCriterion,
        Finder <: NearbySegmentFinder,
        Periods <: Tuple{Vararg{Real}},
    } <: AbstractReconnectionCache
    crit   :: Criterion
    finder :: Finder
    Ls     :: Periods
end

criterion(c::ReconnectionCache) = c.crit
periods(c::ReconnectionCache) = c.Ls

"""
    Reconnections.init_cache(
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament},
        Ls::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) -> AbstractReconnectionCache

Initialise reconnection cache.

Required arguments are a reconnection criterion `crit` and the domain dimensions (or *periods*) `Ls`.

The type of cache will vary depending on the inputs:

- if `crit = NoReconnections()`, then reconnections are disabled and this returns a
  `NullReconnectionCache`;
- otherwise, it returns a `ReconnectionCache`.

In the second case, the detection of reconnection candidates can follow two different strategies:

- if the domain is infinite (default), a naive iteration across all filament segment pairs
is performed;
- if the domain is periodic, an efficient cell lists algorithm is used (see the [`CellLists`](@ref) module).
"""
function init_cache(
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament},
        Ls::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )
    _init_cache(crit, fs, Ls)
end

function _init_cache(crit::ReconnectionCriterion, fs, Ls)
    has_nonperiodic_directions = any(L -> L === Infinity(), Ls)
    r_cut = distance(crit) :: Real
    finder = if has_nonperiodic_directions
        NaiveSegmentFinder(fs, r_cut, Ls)
    else
        CellListSegmentFinder(fs, r_cut, Ls; nsubdiv = Val(1))
    end
    ReconnectionCache(crit, finder, Ls)
end

# Used when reconnections are disabled.
struct NullReconnectionCache <: AbstractReconnectionCache end
_init_cache(::NoReconnections, args...) = NullReconnectionCache()
criterion(::NullReconnectionCache) = NoReconnections()
