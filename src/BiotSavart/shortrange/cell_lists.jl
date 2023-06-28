export CellListsBackend

using StaticArrays: StaticArrays, similar_type, Size
using Static: StaticInt, static, dynamic

"""
    SegmentCellList
    SegmentCellList(
        ::Type{<:AbstractFilament}, rs_cut::NTuple{N, Real}, periods::NTuple{N, Real},
        [nsubdiv::StaticInt = static(1)],
    )

Construct a cell list for dealing with the interaction between filament segments.

Above, `N` is the number of dimensions (usually `N = 3`).

The applied cutoff radius ``r_{\\text{cut}}`` should be much larger than the maximum segment
length ``ℓ``, or should at least account for ``ℓ``.
Basically, if one wants an actual cut-off radius ``r₀``, then the applied cutoff radius passed
to the constructor should be ``r_{\\text{cut}} = r₀ + ℓ``.
Otherwise, a small amount of interactions within ``[r₀ - ℓ, r₀]`` may be missed.

The cutoff radii `rs_cut` (which can be different in each direction) don't need to exactly
divide the domain period `L` into equal pieces, but it's recommended that it does so for
performance reasons.

One can optionally specify a subdivision of the cells by passing `nsubdiv`, which allows to
reduce the number of spurious pairs (i.e. beyond the cutoff radius) as described
[here](https://en.wikipedia.org/wiki/Cell_lists#Improvements). In practice, a value of
`static(2)` can improve performance significantly compared to no subdivision (`static(1)`).

Infinite non-periodic domains (in the sense of `period = Infinity()`) are not supported.
"""
struct SegmentCellList{
        N,  # usually N = 3 (=> 3D space)
        S <: Segment,
        SegmentList <: AbstractVector{S},
        CutoffRadii <: NTuple{N, Real},
        SubDiv <: StaticInt,
        Periods <: NTuple{N, Real},
    }
    segments :: Array{SegmentList, N}  # segments[i, j, k] contains all filament segments inside cell (i, j, k)
    rs_cut   :: CutoffRadii            # cutoff radii (can be different in each direction)
    nsubdiv  :: SubDiv
    Ls       :: Periods
end

function SegmentCellList(
        ::Type{Filament},
        rs_cut_in::NTuple{N, Real},
        Ls::NTuple{N, Real},
        nsubdiv::StaticInt = static(1),
    ) where {N, Filament <: AbstractFilament}
    any(L -> L === Infinity(), Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by CellListsBackend"
    ))

    ndiv = dynamic(nsubdiv)
    rs_cut = map(r -> r / ndiv, rs_cut_in)

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

    SegmentCellList(segs, rs_cut, nsubdiv, Ls)
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
struct CellListsBackend{SubDiv <: StaticInt} <: ShortRangeBackend
    nsubdiv :: SubDiv
end

CellListsBackend() = CellListsBackend(static(1))
CellListsBackend(n::Int) = CellListsBackend(static(n))

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
    (; backend,) = params
    (; rcut,) = params
    (; Ls,) = pc
    # Increase cut-off radius along each direction so that it exactly divides the domain size.
    rs_cut = map(Ls) do L
        L / floor(L / rcut)
    end
    @assert all(≥(rcut), rs_cut)
    cl = SegmentCellList(eltype(fs), rs_cut, Ls, backend.nsubdiv)
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

# Periodic wrapping of array indices
function _wrap_index(inds::AbstractRange, i::Int)
    while i < first(inds)
        i += length(inds)
    end
    while i > last(inds)
        i -= length(inds)
    end
    i
end

_wrap_index(A::AbstractArray{T, N} where {T}, I::NTuple{N}) where {N} =
    map(_wrap_index, axes(A), I)

_wrap_index(A::AbstractArray, I::CartesianIndex) = _wrap_index(A, Tuple(I))

# As of Julia 1.9.1, the @inline is needed to avoid poor performance.
# This seems to be related to the use of Iterators.product (type of `cell_indices`)...
@inline function Base.iterate(it::CellListSegmentIterator, state = nothing)
    (; cl, cell_indices,) = it
    (; segments,) = cl

    if state === nothing  # initial iteration
        cell_index, cell_indices_state = iterate(cell_indices)
        I = _wrap_index(segments, cell_index)
        @inbounds segments_in_current_cell = segments[I...]
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
        I = _wrap_index(segments, cell_index)
        @inbounds segments_in_current_cell = segments[I...]
        ret_segment = iterate(segments_in_current_cell)
        cell_indices_state, segments_in_current_cell, ret_segment
    end

    current_segment, segments_state = ret_segment
    state_next = (cell_indices_state, segments_in_current_cell, segments_state,)
    current_segment, state_next
end

function nearby_segments(c::CellListsCache, x⃗::Vec3)
    (; params, cl,) = c
    (; segments, rs_cut, nsubdiv,) = cl
    (; common,) = params
    (; Ls,) = common
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(segments))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    δ = dynamic(nsubdiv)
    offsets = ntuple(n -> -δ - 1 + n, Val(2δ + 1))  # index offset along each dimension
    inds_per_direction = map(Tuple(I₀)) do i
        map(off -> i + off, offsets)
    end
    cell_indices = Iterators.product(inds_per_direction...)  # iterate over the M³ cells around I₀
    CellListSegmentIterator(cl, cell_indices)
end
