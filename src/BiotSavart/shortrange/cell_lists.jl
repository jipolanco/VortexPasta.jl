export CellListsBackend

using StaticArrays: StaticArrays, similar_type, Size

"""
    SegmentCellList
    SegmentCellList(::Type{<:AbstractFilament}, rcut::Real, periods::NTuple{3, Real})

Construct a cell list for dealing with the interaction between filament segments.

The applied cutoff radius ``r_{\\text{cut}}`` should be much larger than the maximum segment
length ``ℓ``, or should at least account for ``ℓ``.
Basically, if one wants an actual cut-off radius ``r₀``, then the applied cutoff radius passed
to the constructor should be ``r_{\\text{cut}} = r₀ + ℓ``.
Otherwise, a small amount of interactions within ``[r₀ - ℓ, r₀]`` may be missed.

The cutoff radius `rcut` doesn't need to exactly divide the domain period `L` into equal pieces.

Note that infinite non-periodic domains (in the sense of `period = Infinity()`) are not supported.
"""
struct SegmentCellList{
        N,  # usually N = 3 (=> 3D space)
        S <: Segment,
        SegmentList <: AbstractVector{S},
        CutoffRadii <: NTuple{N, Real},
        Periods <: NTuple{N, Real},
    }
    segments :: Array{SegmentList, N}  # segments[i, j, k] contains all filament segments inside cell (i, j, k)
    rs_cut   :: CutoffRadii            # cutoff radii (can be different in each direction)
    Ls       :: Periods
end

function SegmentCellList(
        ::Type{Filament},
        rs_cut::NTuple{N, Real},
        Ls::NTuple{N, Real},
    ) where {N, Filament <: AbstractFilament}
    any(L -> L === Infinity(), Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by CellListsBackend"
    ))

    # Number of cells in each direction.
    # Using `floor` below means that, if `rcut` doesn't exactly divide the domain size L in
    # a given direction, then the *last* cell in that direction will be larger than `rcut`.
    ncells = map(rs_cut, Ls) do rcut, L
        floor(Int, L / rcut)
    end

    S = Segment{Filament}
    @assert isconcretetype(S)

    vempty = Vector{S}(undef, 0)
    # vempty = StructVector{S}(undef, 0)  # not necessarily faster than a regular Vector in this case
    SegmentList = typeof(vempty)
    segs = Array{SegmentList, N}(undef, ncells)

    for i ∈ eachindex(segs)
        segs[i] = copy(vempty)
    end

    SegmentCellList(segs, rs_cut, Ls)
end

function Base.empty!(cl::SegmentCellList)
    for v ∈ cl.segments
        empty!(v)
    end
    cl
end

@inline function determine_cell_index(x, rcut, L, N)
    while x < 0
        x += L
    end
    while x ≥ L
        x -= L
    end
    clamp(1 + floor(Int, x / rcut), 1, N)  # make sure the index is in 1:N
end

function add_segment!(cl::SegmentCellList{N, S}, seg::S) where {N, S <: Segment}
    (; segments, rs_cut, Ls,) = cl
    x⃗ = Filaments.midpoint(seg)
    inds = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(segments))
    I = CartesianIndex(inds)
    @inbounds push!(segments[I], seg)
    cl
end

function assign_cells!(cl::SegmentCellList, f::AbstractFilament)
    for s ∈ segments(f)
        add_segment!(cl, s)
    end
    cl
end

function assign_cells!(cl::SegmentCellList, fs::AbstractVector{<:AbstractFilament})
    empty!(cl)
    for f ∈ fs
        assign_cells!(cl, f)
    end
    cl
end

# ================================================================================ #

"""
    CellListsBackend <: ShortRangeBackend

Compute short-range interactions using the cell lists algorithm.

This backend can be significantly faster than the [`NaiveShortRangeBackend`](@ref) when the
cutoff radius `rcut` is much smaller than the domain period `L` (roughly when `rcut ≲ L / 10`).

Future improvements may further increase the performance of this backend.

This backend does not support non-periodic domains.

See [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for details.
"""
struct CellListsBackend <: ShortRangeBackend end

struct CellListsCache{
        CellList <: SegmentCellList,
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
    (; rcut,) = params
    (; Ls,) = pc
    # Increase cut-off radius along each direction so that it exactly divides the domain size.
    rs_cut = map(Ls) do L
        L / floor(L / rcut)
    end
    @assert all(≥(rcut), rs_cut)
    cl = SegmentCellList(eltype(fs), rs_cut, Ls)
    CellListsCache(cl, params, to)
end

function set_filaments!(c::CellListsCache, fs)
    assign_cells!(c.cl, fs)
    c
end

struct CellListSegmentIterator{
        S <: Segment,
        N,
        CellList <: SegmentCellList{N, S},
        CellIndices,
    } <: NearbySegmentIterator{S}
    cl           :: CellList
    cell_indices :: CellIndices  # iterator over indices of cells to be visited
    function CellListSegmentIterator(cl::SegmentCellList{N, S}, inds) where {N, S}
        new{S, N, typeof(cl), typeof(inds)}(cl, inds)
    end
end

# First iteration.
# As of Julia 1.9.1, the @inline is needed to avoid poor performance, which seems to be
# related to the use of Iterators.product...
@inline function Base.iterate(it::CellListSegmentIterator)
    (; cl, cell_indices,) = it
    (; segments,) = cl

    cell_index, cell_indices_state = iterate(cell_indices)
    segments_in_current_cell = segments[cell_index...] :: AbstractVector{<:Segment}
    ret_segment = iterate(segments_in_current_cell)

    # If the initial cell has no segments, continue iterating until we find a cell with
    # segments, or until we've visited all cells.
    while ret_segment === nothing  # case of a cell with no segments
        ret = _jump_to_next_cell(segments, cell_indices, cell_indices_state)
        ret === nothing && return nothing  # we're done iterating over cells
        cell_indices_state, segments_in_current_cell, ret_segment = ret
    end

    current_segment, segments_state = ret_segment
    state_next = (cell_indices_state, segments_in_current_cell, segments_state,)
    current_segment, state_next
end

function Base.iterate(it::CellListSegmentIterator, state)
    (; cl, cell_indices,) = it
    (; segments,) = cl
    (cell_indices_state, segments_in_current_cell, segments_state,) = state

    # 1. Try to keep iterating over the segments of the current cell.
    ret_segment = iterate(segments_in_current_cell, segments_state)
    if ret_segment !== nothing
        current_segment, segments_state = ret_segment
        state_next = (cell_indices_state, segments_in_current_cell, segments_state,)
        return current_segment, state_next
    end

    # 2. We're done iterating over the current cell, so we jump to the next non-empty cell.
    while ret_segment === nothing
        ret = _jump_to_next_cell(segments, cell_indices, cell_indices_state)
        ret === nothing && return nothing  # we're done iterating over cells
        cell_indices_state, segments_in_current_cell, ret_segment = ret
    end

    current_segment, segments_state = ret_segment
    state_next = (cell_indices_state, segments_in_current_cell, segments_state,)
    current_segment, state_next
end

function _jump_to_next_cell(segments, cell_indices, cell_indices_state)
    ret_cell = iterate(cell_indices, cell_indices_state)
    ret_cell === nothing && return nothing  # we're done iterating over cells
    cell_index, cell_indices_state = ret_cell
    segments_in_current_cell = segments[cell_index...]
    ret_segment = iterate(segments_in_current_cell)
    cell_indices_state, segments_in_current_cell, ret_segment
end

function nearby_segments(c::CellListsCache, x⃗::Vec3)
    (; params, cl,) = c
    (; segments, rs_cut,) = cl
    (; common,) = params
    (; Ls,) = common
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(segments))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    offsets = (-1, 0, 1)  # index offset along each dimension
    inds_per_direction = map(Tuple(I₀), size(segments)) do i, N
        map(offsets) do δi
            j = i + δi
            if j ≤ 0
                j += N  # periodic wrapping
            elseif j > N
                j -= N
            end
            j  # absolute index of cell along one dimension
        end
    end
    cell_indices = Iterators.product(inds_per_direction...)  # iterate over the M³ cells around I₀
    CellListSegmentIterator(cl, cell_indices)
end
