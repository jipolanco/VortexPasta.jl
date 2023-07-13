export CellListsBackend

using ..CellLists: CellLists, PeriodicCellList, add_element!, static,
                   determine_cell_index  # TODO remove?

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
        add_element!(cl, s)
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

struct CellListSegmentIterator{
        S <: Segment,
        N,
        CellList <: PeriodicCellList{N, S},
        CellIndices,
    } <: NearbySegmentIterator{S}
    cl           :: CellList
    cell_indices :: CellIndices  # iterator over indices of cells to be visited
    function CellListSegmentIterator(cl::PeriodicCellList{N, S}, inds) where {N, S}
        new{S, N, typeof(cl), typeof(inds)}(cl, inds)
    end
end

# As of Julia 1.9.1, the @inline is needed to avoid poor performance and spurious allocations.
@inline function Base.iterate(it::CellListSegmentIterator, state = nothing)
    (; cl, cell_indices,) = it
    (; data,) = cl

    if state === nothing  # initial iteration
        cell_index, cell_indices_state = iterate(cell_indices)
        @inbounds segments_in_current_cell = data[cell_index]
        ret_segment = iterate(segments_in_current_cell)  # get first segment of first cell (or `nothing`, if the cell is empty)
    else
        (cell_indices_state, segments_in_current_cell, segments_state,) = state
        ret_segment = iterate(segments_in_current_cell, segments_state)  # advance to next segment (or `nothing`, if we're done with this cell)
    end

    # 1. Try to keep iterating over the segments of the current cell.
    if ret_segment !== nothing
        current_segment, segments_state = ret_segment
        state_next = (cell_indices_state, segments_in_current_cell, segments_state,)
        return current_segment, state_next
    end

    # 2. We're done iterating over the current cell, so we jump to the next non-empty cell.
    while ret_segment === nothing
        ret_cell = iterate(cell_indices, cell_indices_state)
        ret_cell === nothing && return nothing  # we're done iterating over cells
        cell_index, cell_indices_state = ret_cell
        @inbounds segments_in_current_cell = data[cell_index]
        ret_segment = iterate(segments_in_current_cell)
        cell_indices_state, segments_in_current_cell, ret_segment
    end

    current_segment, segments_state = ret_segment
    state_next = (cell_indices_state, segments_in_current_cell, segments_state,)
    current_segment, state_next
end

function nearby_segments(c::CellListsCache, x⃗::Vec3)
    (; params, cl,) = c
    (; data, rs_cut,) = cl
    (; common,) = params
    (; Ls,) = common
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(data))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    M = CellLists.subdivisions(cl)
    cell_indices = CartesianIndices(
        map(i -> (i - M):(i + M), Tuple(I₀))
    )
    CellListSegmentIterator(cl, cell_indices)
end
