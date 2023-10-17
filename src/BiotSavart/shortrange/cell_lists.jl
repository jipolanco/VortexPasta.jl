export CellListsBackend

using ..FindNearbySegments: CellListSegmentFinder
using ..CellLists: PeriodicCellList  # for docs only

"""
    CellListsBackend <: ShortRangeBackend
    CellListsBackend(nsubdiv::Int = 1)

Compute short-range interactions using the cell lists algorithm.

This backend can be significantly faster than the [`NaiveShortRangeBackend`](@ref) when the
cutoff radius `rcut` is much smaller than the domain period `L` (roughly when `rcut ≲ L / 10`).

Optionally, one can choose to subdivide each cell (of size `≈ rcut`) onto `nsubdiv`
subcells. In practice, a value of `2` or `3` can significantly improve performance compared
to no subdivision (`1`).

This backend does not support non-periodic domains.

See [`PeriodicCellList`](@ref) and [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for
more details.
"""
struct CellListsBackend{M} <: ShortRangeBackend end
CellListsBackend(n::Int = 1) = CellListsBackend{n}()

subdivisions(::CellListsBackend{M}) where {M} = M

struct CellListsCache{
        Finder <: CellListSegmentFinder,
        Params <: ParamsShortRange{<:CellListsBackend},
        Timer <: TimerOutput,
    } <: ShortRangeCache
    finder :: Finder
    params :: Params
    to     :: Timer
end

function init_cache_short(
        pc::ParamsCommon, params::ParamsShortRange{<:CellListsBackend},
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    )
    (; backend, rcut,) = params
    (; Ls,) = pc
    nsubdiv = Val(subdivisions(backend))
    finder = CellListSegmentFinder(fs, rcut, Ls; nsubdiv)
    CellListsCache(finder, params, to)
end
