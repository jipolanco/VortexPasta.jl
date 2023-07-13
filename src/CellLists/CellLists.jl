"""
    CellLists

Module implementing the cell lists algorithm over N-dimensional periodic domains.
"""
module CellLists

export PeriodicCellList, CellListIterator, static, nearby_elements

using ..PaddedArrays: PaddedArray, pad_periodic!
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

Base.size(cl::PeriodicCellList) = size(cl.data)
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

"""
    Base.empty!(cl::PeriodicCellList) -> cl

Remove all elements from the cell list.
"""
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

"""
    add_element!(cl::PeriodicCellList{N, T}, el::T)

Add element to the cell list.

Determines the cell associated to the element and then appends the element to that cell.
"""
function add_element!(cl::PeriodicCellList{N, T}, el::T) where {N, T}
    (; data, rs_cut, Ls, to_coordinate,) = cl
    x⃗ = to_coordinate(el)
    inds = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(cl))
    I = CartesianIndex(inds)
    @inbounds cell = data[I]
    push!(cell, el)  # this can allocate if we don't put a `sizehint!` somewhere (which we do!)
    cl
end

struct CellListIterator{
        T, N,
        CellList <: PeriodicCellList{N, T},
        CellIndices,
    }
    cl           :: CellList
    cell_indices :: CellIndices  # iterator over indices of cells to be visited
    function CellListIterator(cl::PeriodicCellList{N, T}, inds) where {N, T}
        new{T, N, typeof(cl), typeof(inds)}(cl, inds)
    end
end

# These are needed e.g. by collect(it::CellListIterator)
Base.IteratorSize(::Type{<:CellListIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:CellListIterator{T}}) where {T} = T

# As of Julia 1.9.1, the @inline is needed to avoid poor performance and spurious allocations.
@inline function Base.iterate(it::CellListIterator, state = nothing)
    (; cl, cell_indices,) = it
    (; data,) = cl

    if state === nothing  # initial iteration
        cell_index, cell_indices_state = iterate(cell_indices)
        @inbounds elements_in_current_cell = data[cell_index]
        ret_element = iterate(elements_in_current_cell)  # get first element of first cell (or `nothing`, if the cell is empty)
    else
        (cell_indices_state, elements_in_current_cell, elements_state,) = state
        ret_element = iterate(elements_in_current_cell, elements_state)  # advance to next element (or `nothing`, if we're done with this cell)
    end

    # 1. Try to keep iterating over the elements of the current cell.
    if ret_element !== nothing
        current_element, elements_state = ret_element
        state_next = (cell_indices_state, elements_in_current_cell, elements_state,)
        return current_element, state_next
    end

    # 2. We're done iterating over the current cell, so we jump to the next non-empty cell.
    while ret_element === nothing
        ret_cell = iterate(cell_indices, cell_indices_state)
        ret_cell === nothing && return nothing  # we're done iterating over cells
        cell_index, cell_indices_state = ret_cell
        @inbounds elements_in_current_cell = data[cell_index]
        ret_element = iterate(elements_in_current_cell)
        cell_indices_state, elements_in_current_cell, ret_element
    end

    current_element, elements_state = ret_element
    state_next = (cell_indices_state, elements_in_current_cell, elements_state,)
    current_element, state_next
end

"""
    nearby_elements(cl::PeriodicCellList{N}, x⃗) -> CellListIterator

Return an iterator over the elements that are sufficiently close to the point `x⃗`.

The iterator returns the elements which are likely to be within the cutoff radius of the
point `x⃗`. More precisely, it returns elements in the same cell as `x⃗` as well as in
neighbouring cells.

Here x⃗ should be a coordinate, usually represented by an `SVector{N}` or an `NTuple{N}`.
"""
function nearby_elements(cl::PeriodicCellList{N}, x⃗) where {N}
    length(x⃗) == N || throw(DimensionMismatch(lazy"wrong length of coordinate: x⃗ = $x⃗ (expected length is $N)"))
    (; rs_cut, Ls,) = cl
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(cl))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    M = subdivisions(cl)
    cell_indices = CartesianIndices(
        map(i -> (i - M):(i + M), Tuple(I₀))
    )
    CellListIterator(cl, cell_indices)
end

end
