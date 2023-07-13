export CellListsBackend

using ..CellLists:
    CellLists,
    PeriodicCellList,
    CellListIterator,
    static

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
        CellList <: PeriodicCellList,
        Params <: ParamsShortRange,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    cl     :: CellList
    params :: Params
    to     :: Timer
end

function init_cache_short(
        pc::ParamsCommon, params::ParamsShortRange{<:CellListsBackend},
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    )
    (; backend,) = params
    (; rcut,) = params
    (; Ls,) = pc
    # Increase cut-off radius along each direction so that it exactly divides the domain size.
    rs_cut = map(Ls) do L
        L / floor(L / rcut)
    end
    @assert all(≥(rcut), rs_cut)
    M = subdivisions(backend)
    Filament = eltype(fs)
    S = Segment{Filament}
    @assert isconcretetype(S)
    # Construct cell list of filament segments.
    # We use the `midpoint` function to associate a coordinate to each segment.
    cl = PeriodicCellList(S, rs_cut, Ls, static(M); to_coordinate = Filaments.midpoint)
    CellListsCache(cl, params, to)
end

function assign_cells!(cl::PeriodicCellList, f::AbstractFilament)
    for s ∈ segments(f)
        CellLists.add_element!(cl, s)
    end
    cl
end

function assign_cells!(cl::PeriodicCellList, fs::AbstractVector{<:AbstractFilament})
    empty!(cl)
    for f ∈ fs
        assign_cells!(cl, f)
    end
    cl
end

function set_filaments!(c::CellListsCache, fs)
    assign_cells!(c.cl, fs)
    c
end

nearby_segments(c::CellListsCache, x⃗::Vec3) = CellLists.nearby_elements(c.cl, x⃗)
