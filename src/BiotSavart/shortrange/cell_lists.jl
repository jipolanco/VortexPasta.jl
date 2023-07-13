export CellListsBackend

using ..BasicTypes: PaddedArray, pad_periodic!
using StaticArrays: StaticArrays, similar_type, Size
using Static: StaticInt, static, dynamic

"""
    PeriodicCellList{N, T}
    PeriodicCellList(
        ::Type{T}, rs_cut::NTuple{N, Real}, periods::NTuple{N, Real},
        [nsubdiv::StaticInt = static(1)];
        [to_coordinate::Function = identity],
    )

Construct a cell list for dealing with pair interactions.

Above, `N` is the number of spatial dimensions, and `T` is the type of each element.
In the simplest cases, `T` can simply describe a coordinate in N-dimensional space
(e.g. `T = SVector{N, Float64}`). One can also deal with more complicated elements which
include more information. As an example, see further below for how to deal with filament
segments in 3D space.

The cutoff radii `rs_cut` (which can be different in each direction) don't need to exactly
divide the domain period `L` into equal pieces, but it's recommended that it does so for
performance reasons.

Optionally, one can choose to subdivide each cell (of size `≈ rcut`) onto `nsubdiv`
subcells. This can significantly improve performance, since it allows to discard some
spurious pair interactions (i.e. beyond the chosen cutoff radius) as described
[here](https://en.wikipedia.org/wiki/Cell_lists#Improvements). In practice, a value of
`2` or `3` can significantly improve performance compared to no subdivision (`1`).

Infinite non-periodic domains (in the sense of `period = Infinity()`) are not supported.

# Dealing with filament segments

One of the possible uses of `PeriodicCellList` is to classify filament segments (which are
typically shorter than the cutoff radius) according to their spatial location. In that case,
`T` is not a simple coordinate, but may contain more information including things like (1)
the filament the segment belongs to, and (2) the location of the segment within the filament.
As there is no unique way of associating a coordinate to a segment, one should pass the
`to_coordinate` argument which "converts" the segment to a coordinate in space. For instance,
the passed `to_coordinate` function may return the midpoint of the segment, which will be
used to determine the cell associated to the segment.

The applied cutoff radius ``r_{\\text{cut}}`` should be much larger than the maximum segment
length ``ℓ``, or should at least account for ``ℓ``.
Basically, if one wants an actual cut-off radius ``r₀``, then the applied cutoff radius passed
to the constructor should be ``r_{\\text{cut}} = r₀ + ℓ``.
Otherwise, a small amount of interactions within ``[r₀ - ℓ, r₀]`` may be missed.
"""
struct PeriodicCellList{
        N,  # spatial dimension
        T,  # type of "element": can be a coordinate in 3D space, or some other element which includes some coordinate information (e.g. a filament segment)
        ElementList <: AbstractVector{T},
        M,  # number of cell subdivisions (≥ 1)
        CellList <: PaddedArray{M, ElementList, N},  # array of cells, each containing a list of elements
        ToCoordFunc <: Function,
        CutoffRadii <: NTuple{N, Real},
        Periods <: NTuple{N, Real},
    }
    data          :: CellList     # data[i, j, k] contains all "elements" inside cell (i, j, k)
    to_coordinate :: ToCoordFunc  # convert "element" to N-D coordinate
    rs_cut        :: CutoffRadii  # cutoff radii (can be different in each direction)
    nsubdiv       :: StaticInt{M}
    Ls            :: Periods
end

subdivisions(::PeriodicCellList{A, B, C, M}) where {A, B, C, M} = M

function PeriodicCellList(
        ::Type{T},
        rs_cut_in::NTuple{N, Real},
        Ls::NTuple{N, Real},
        nsubdiv::StaticInt = static(1);
        to_coordinate::F = identity,
    ) where {N, T, F <: Function}
    any(isinf, Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by PeriodicCellList"
    ))

    M = dynamic(nsubdiv)
    rs_cut = map(r -> r / M, rs_cut_in)

    # Number of cells in each direction.
    # Using `floor` below means that, if `rcut` doesn't exactly divide the domain size L in
    # a given direction, then the *last* cell in that direction will be larger than `rcut`.
    ncells = map(rs_cut, Ls) do rcut, L
        floor(Int, L / rcut)
    end

    all(≥(2M), ncells) || error(
        lazy"""number of cells $ncells is too small for periodic padding.
               Try reducing the cutoff radius (got rs_cut = $rs_cut_in)."""
    )

    @assert isconcretetype(T)
    vempty = Vector{T}(undef, 0)
    # vempty = StructVector{S}(undef, 0)  # not necessarily faster than a regular Vector in this case
    ElementList = typeof(vempty)
    data_dims = ncells .+ 2M  # add 2M ghost cells in each direction
    data_raw = Array{ElementList, N}(undef, data_dims)

    data = PaddedArray{M}(data_raw)
    @assert size(data) == ncells

    # Initialise cells inside the domain (i.e. not including ghost cells)
    for I ∈ CartesianIndices(data)
        data[I] = copy(vempty)
    end

    # Pad array periodically. Note that this needs to be done only once (and not whenever
    # filaments are added), since we're copying array references ("pointers"), so modifying
    # a "central" cell will also modify its corresponding ghost cell if it has one.
    pad_periodic!(data)

    PeriodicCellList(data, to_coordinate, rs_cut, nsubdiv, Ls)
end

function Base.empty!(cl::PeriodicCellList)
    for v ∈ cl.data
        n = length(v)
        empty!(v)
        # Heuristic to reduce allocations in `push!` when refilling the cells.
        sizehint!(v, n < 8 ? 8 : nextpow(2, n))
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

function add_element!(cl::PeriodicCellList{N, T}, el::T) where {N, T}
    (; data, rs_cut, Ls, to_coordinate,) = cl
    x⃗ = to_coordinate(el)
    inds = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(data))
    I = CartesianIndex(inds)
    @inbounds cell = data[I]
    push!(cell, el)  # this can allocate if we don't put a `sizehint!` somewhere (which we do!)
    cl
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

# ================================================================================ #

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
    M = subdivisions(cl)
    cell_indices = CartesianIndices(
        map(i -> (i - M):(i + M), Tuple(I₀))
    )
    CellListSegmentIterator(cl, cell_indices)
end
