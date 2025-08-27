"""
    CellLists

Module implementing the cell lists algorithm over ``N``-dimensional periodic domains.

See the [Wikipedia article](https://en.wikipedia.org/wiki/Cell_lists) for some details.
"""
module CellLists

# The implementation is based on
# https://aiichironakano.github.io/cs596/01-1LinkedListCell.pdf with a few differences:
# - the treatment of periodicity is different (we use ghost cells)
# - they describe the case with nsubdiv = 1, while we allow nsubdiv ≥ 1

export PeriodicCellList, static, nearby_elements

using ..PaddedArrays: PaddedArray, pad_periodic!
using Static: StaticInt, static, dynamic

# This is used either to mean that a cell has no elements, or that an element is the last
# element within a given cell. This assumes one-based indexing, so that 0 is an invalid index!
const EMPTY = 0

## ================================================================================ ##
## Iteration over a single cell

struct CellIterator{
        T,
        ElementList <: AbstractVector{T},
        IndexType <: Integer
    }
    elements   :: ElementList        # contains all elements in all cells
    head_index :: IndexType          # index of first element in cell
    next_index :: Vector{IndexType}  # allows to get the rest of the elements in cell
end

Base.IteratorSize(::Type{<:CellIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:CellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:CellIterator{T}}) where {T} = T

function Base.iterate(it::CellIterator, n = it.head_index)
    (; elements, next_index,) = it
    n == EMPTY && return nothing  # no more elements in this cell
    el = @inbounds elements[n]
    @inbounds el, next_index[n]
end

## ================================================================================ ##
## Cell list type definition

"""
    PeriodicCellList{N, T}
    PeriodicCellList(
        ::Type{T}, rs_cut::NTuple{N, Real}, periods::NTuple{N, Real},
        [nsubdiv = static(1)];
        [to_coordinate::Function = identity],
    )

Construct a cell list for dealing with pair interactions.

Above, `N` is the number of spatial dimensions, and `T` is the type of each element.
In the simplest cases, `T` can simply describe a coordinate in ``N``-dimensional space
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
For convenience, it can be passed as `nsubdiv = Val(M)` or as `nsubdiv = static(M)`.

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
        HeadArray <: PaddedArray{M, Int, N},  # array of cells, each containing a list of elements
        RefBool <: Ref{Bool},
        ToCoordFunc <: Function,
        CutoffRadii <: NTuple{N, Real},
        Periods <: NTuple{N, Real},
    }
    elements      :: ElementList  # contains all "elements" (unsorted) [length Np]
    head_indices  :: HeadArray    # array pointing to the index of the first element in each cell (or EMPTY if no elements in that cell)
    next_index    :: Vector{Int}  # points from one element to the next element (or EMPTY if this is the last element) [length Np]
    isready       :: RefBool
    to_coordinate :: ToCoordFunc  # convert "element" to N-D coordinate
    rs_cut        :: CutoffRadii  # cutoff radii (can be different in each direction)
    nsubdiv       :: StaticInt{M}
    Ls            :: Periods
end

Base.size(cl::PeriodicCellList) = size(cl.head_indices)
subdivisions(::PeriodicCellList{A, B, C, M}) where {A, B, C, M} = M

@inline PeriodicCellList(::Type{T}, rs, Ls, nsubdiv::Val{M}; kws...) where {T, M} =
    PeriodicCellList(T, rs, Ls, static(M); kws...)

# Returns maximum possible cut-off distance r_cut for M subdivisions and a domain of period L.
max_cutoff_distance(M::Integer, L::AbstractFloat) = oftype(L, M / (2M + 1)) * L

function PeriodicCellList(
        ::Type{T},
        rs_cut::NTuple{N, Real},
        Ls::NTuple{N, Real},
        nsubdiv::StaticInt = static(1);
        to_coordinate::F = identity,
    ) where {N, T, F <: Function}
    any(isinf, Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by PeriodicCellList"
    ))

    M = dynamic(nsubdiv)

    # Determine cell sizes.
    rs_cell = map(rs_cut, Ls) do r, L
        r′ = r / M
        # Increase cell size so that it exactly divides the domain size.
        L / floor(L / r′)
    end

    # Number of cells in each direction.
    # Using `floor` below means that, if `r_cell` doesn't exactly divide the domain size L in
    # a given direction, then the *last* cell in that direction will be larger than `r_cell`.
    ncells = map(rs_cell, Ls) do r_cell, L
        floor(Int, L / r_cell)
    end

    # When M = 1, the number of cells in each direction should be at least 3 to avoid
    # repeating pair interactions (due to periodicity).
    # More generally, for any M, the number of cells should be at least 2M + 1.
    # In the end, this corresponds to the condition r_cut_in/L < M / (2M + 1).
    all(≥(2M + 1), ncells) || error(
        lazy"""cell lists: the cut-off distance r_cut is too large for periodic padding.
               It should satisfy r_cut/L ≤ nsubdiv / (2 * nsubdiv + 1) = $(M / (2M + 1));
               got rs_cut/Ls = $(rs_cut ./ Ls)."""
    )

    elements = Vector{T}(undef, 0)
    # elements = StructVector{T}(undef, 0)  # not necessarily faster than a regular Vector in this case

    IndexType = Int
    head_dims = ncells .+ 2M  # add 2M ghost cells in each direction
    head_raw = Array{IndexType}(undef, head_dims)
    fill!(head_raw, EMPTY)

    head_indices = PaddedArray{M}(head_raw)
    @assert size(head_indices) == ncells

    next_index = Vector{Int}(undef, 0)
    isready = Ref(true)  # technically we're ready for iterating (with zero elements)

    Base.require_one_based_indexing(elements)  # assumed since EMPTY = 0

    PeriodicCellList(
        elements, head_indices, next_index, isready, to_coordinate, rs_cell, nsubdiv, Ls,
    )
end

"""
    Base.empty!(cl::PeriodicCellList) -> cl

Remove all elements from the cell list.
"""
function Base.empty!(cl::PeriodicCellList)
    empty!(cl.elements)
    empty!(cl.next_index)
    fill!(parent(cl.head_indices), EMPTY)  # this can be slow with too many cells!
    cl
end

function Base.sizehint!(cl::PeriodicCellList, Np)
    sizehint!(cl.elements, Np)
    sizehint!(cl.next_index, Np)
    cl
end

@inline function determine_cell_index(x, rcut, L, N)
    while x < 0
        x += L
    end
    while x ≥ L
        x -= L
    end
    # Here unsafe_trunc(Int, ⋅) is used instead of floor(Int, ⋅) because it should be
    # faster.
    # For non-negative values, both should give the same result.
    # The unsafe_trunc function generally uses a single intrinsic CPU instruction and never
    # throws errors. It can silently give a wrong result if the values are not representable
    # by an Int, but that will never be the case in practice here (since 0 ≤ x/rcut < L/rcut
    # and L/rcut is very small compared to typemax(Int) = 2^63 - 1).
    clamp(1 + unsafe_trunc(Int, x / rcut), 1, N)  # make sure the index is in 1:N
end

"""
    add_element!(cl::PeriodicCellList{N, T}, el::T, [x⃗])

Add element to the cell list.

Determines the cell associated to the element and then appends the element to that cell.

Optionally, one may pass the coordinate location ``x⃗`` associated to the element.
Otherwise, it will be obtained from the element according to

    x⃗ = to_coordinate(el)

where `to_coordinate` corresponds to the keyword argument of [`PeriodicCellList`](@ref).
"""
function add_element!(cl::PeriodicCellList{N, T}, el::T) where {N, T}
    x⃗ = cl.to_coordinate(el)
    add_element!(cl, el, x⃗)
end

function add_element!(cl::PeriodicCellList{N, T}, el::T, x⃗) where {N, T}
    (; rs_cut, Ls, elements, next_index, head_indices,) = cl
    cl.isready[] = false      # we're not ready for iterating
    @assert length(elements) == length(next_index)
    push!(elements, el)       # add this element to element vector
    inds = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(cl))
    I = CartesianIndex(inds)
    n = lastindex(elements)                        # index of the new element
    push!(next_index, @inbounds(head_indices[I]))  # the old head now comes after the new element
    @inbounds head_indices[I] = n                  # the new element is the new head
    cl
end

"""
    CellLists.finalise_cells!(cl::PeriodicCellList)

"Finalise" cells before iterating over its elements.

This function performs operations needed to iterate over elements of the cell lists, such as
filling ghost cells.
It must be called once after all elements have been added to the cell list using
[`add_element!`](@ref) and before one starts iterating using [`nearby_elements`](@ref).
"""
function finalise_cells!(cl::PeriodicCellList)
    cl.isready[] && return cl
    pad_periodic!(cl.head_indices)  # fill ghost cells for periodicity (can be slow...)
    cl.isready[] = true
    cl
end

"""
    CellLists.set_elements!(get_element::Function, cl::PeriodicCellList, Np::Integer)

Set all elements of the cell list.

Here `get_element(n::Integer)` is a function that returns a single element given its index `n`.
It must take and index `n` in `1:Np` (this assumes one-based indexing!).
`Np` is the total number of elements.

This function resets the cell list, removing all previously existent points.
It can be used as a replacement for [`empty(::PeriodicCellList)`](@ref) + [`add_element!`](@ref) + [`finalise_cells!`](@ref).
"""
function set_elements!(get_element::F, cl::PeriodicCellList{N, T}, Np::Integer) where {F <: Function, N, T}
    (; elements, next_index, head_indices, Ls, rs_cut,) = cl
    fill!(parent(head_indices), EMPTY)  # this can be slow with too many cells?
    resize!(elements, Np)
    resize!(next_index, Np)
    Base.require_one_based_indexing(elements)
    Base.require_one_based_indexing(next_index)
    @inbounds for n in 1:Np
        el = @inline get_element(n)
        elements[n] = el
        x⃗ = @inline cl.to_coordinate(el)
        inds = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(cl))
        I = CartesianIndex(inds)
        next_index[n] = head_indices[I]  # the old head now comes after the new element
        head_indices[I] = n              # the new element is the new head
    end
    cl.isready[] = true
    cl
end

## ================================================================================ ##
## Iteration over elements in cell lists

struct MultiCellIterator{
        T, N,
        ElementList <: AbstractVector{T},
        IndexType,
        HeadArray <: AbstractArray{IndexType, N},
        CellIndices <: CartesianIndices{N},
    }
    elements     :: ElementList        # contains all elements in all cells
    head_indices :: HeadArray          # array pointing to the index of the first element in each cell
    next_index   :: Vector{IndexType}  # allows to get the rest of the elements in cell
    cell_indices :: CellIndices        # indices of cells that will be visited
end

Base.IteratorSize(::Type{<:MultiCellIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:MultiCellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:MultiCellIterator{T}}) where {T} = T

@inline function Base.iterate(it::MultiCellIterator)
    icell = firstindex(it.cell_indices)
    I = @inbounds it.cell_indices[icell]
    n = @inbounds it.head_indices[I]
    iterate(it, (icell, n))
end

@inline function Base.iterate(it::MultiCellIterator, state)
    (; elements, head_indices, next_index, cell_indices,) = it
    icell, n = state
    while n == EMPTY
        if icell == lastindex(cell_indices)
            return nothing  # if this was the last cell, stop iterating
        end
        icell += 1  # jump to the next cell
        I = @inbounds cell_indices[icell]
        n = @inbounds head_indices[I]
    end
    el = @inbounds elements[n]
    @inbounds el, (icell, next_index[n])
end

"""
    nearby_elements(cl::PeriodicCellList{N}, x⃗)

Return an iterator over the elements that are sufficiently close to the point `x⃗`.

The iterator returns the elements which are likely to be within the cutoff radius of the
point `x⃗`. More precisely, it returns elements in the same cell as `x⃗` as well as in
neighbouring cells.

Here `x⃗` should be a coordinate, usually represented by an `SVector{N}` or an `NTuple{N}`.
"""
function nearby_elements(cl::PeriodicCellList{N}, x⃗) where {N}
    length(x⃗) == N || throw(DimensionMismatch(lazy"wrong length of coordinate: x⃗ = $x⃗ (expected length is $N)"))
    (; rs_cut, Ls,) = cl
    cl.isready[] || error("one must call `finalise_cells!` on the `PeriodicCellList` before iterating using `nearby_elements`")
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(cl))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    M = subdivisions(cl)
    cell_indices = CartesianIndices(
        map(i -> (i - M):(i + M), Tuple(I₀))
    )
    # iters = (CellIterator(cl.elements, cl.head_indices[I], cl.next_index) for I ∈ cell_indices)
    # it = Iterators.flatten(iters)  # this gives eltype(it) == Any (but that doesn't seem to affect performance...)
    it = MultiCellIterator(cl.elements, cl.head_indices, cl.next_index, cell_indices)  # slightly slower than Iterators.flatten (but with known eltype)
    # eltype(it)
    it
end

end
