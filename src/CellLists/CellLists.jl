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

using ..PaddedArrays: PaddedArrays, PaddedArray, pad_periodic!
using Atomix: Atomix
using Static: StaticInt, static, dynamic

# This is used either to mean that a cell has no elements, or that an element is the last
# element within a given cell. This assumes one-based indexing, so that 0 is an invalid index!
const EMPTY = 0

## ================================================================================ ##
## Cell list type definition

"""
    PeriodicCellList{N}
    PeriodicCellList(
        rs_cut::NTuple{N, Real}, periods::NTuple{N, Real},
        [nsubdiv = Val(1)];
        [to_coordinate::Function = identity],
    )

Construct a cell list for dealing with pair interactions in `N` dimensions.

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
"""
struct PeriodicCellList{
        N,  # spatial dimension
        M,  # number of cell subdivisions (≥ 1)
        HeadArray <: PaddedArray{M, Int, N},  # array of cells, each containing a list of elements
        CellDims <: NTuple{N, Real},
        Periods <: NTuple{N, Real},
    }
    head_indices  :: HeadArray    # array pointing to the index of the first element in each cell (or EMPTY if no elements in that cell)
    next_index    :: Vector{Int}  # points from one element to the next element (or EMPTY if this is the last element) [length Np]
    rs_cell       :: CellDims     # dimensions of a cell (can be different in each direction)
    nsubdiv       :: StaticInt{M}
    Ls            :: Periods
end

Base.size(cl::PeriodicCellList) = size(cl.head_indices)
subdivisions(::PeriodicCellList{N, M}) where {N, M} = M

@inline PeriodicCellList(rs, Ls, nsubdiv::Val{M}; kws...) where {M} = PeriodicCellList(rs, Ls, static(M); kws...)

# Returns maximum possible cut-off distance r_cut for M subdivisions and a domain of period L.
max_cutoff_distance(M::Integer, L::AbstractFloat) = oftype(L, M / (2M + 1)) * L

function PeriodicCellList(
        rs_cut::NTuple{N, Real},
        Ls::NTuple{N, Real},
        nsubdiv::StaticInt = static(1),
    ) where {N}
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

    IndexType = Int
    head_dims = ncells .+ 2M  # add 2M ghost cells in each direction
    head_raw = Array{IndexType}(undef, head_dims)
    fill!(head_raw, EMPTY)

    head_indices = PaddedArray{M}(head_raw)
    @assert size(head_indices) == ncells

    next_index = Vector{IndexType}(undef, 0)

    PeriodicCellList(head_indices, next_index, rs_cell, nsubdiv, Ls)
end

"""
    Base.empty!(cl::PeriodicCellList) -> cl

Remove all elements from the cell list.
"""
function Base.empty!(cl::PeriodicCellList)
    empty!(cl.next_index)
    fill!(parent(cl.head_indices), EMPTY)  # this can be slow with too many cells!
    cl
end

function Base.sizehint!(cl::PeriodicCellList, Np)
    sizehint!(cl.next_index, Np)
    cl
end

# TODO: not sure this is fast on GPUs...
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
    CellLists.set_elements!(cl::PeriodicCellList, xp::AbstractVector)

Set all elements of the cell list.

Here `xp` is a vector of spatial locations.

This function resets the cell list, removing all previously existent points.
"""
function set_elements!(get_coordinate::F, cl::PeriodicCellList, xp::AbstractVector) where {F}
    (; next_index, head_indices, Ls, rs_cell,) = cl
    head_indices_data = parent(head_indices)   # full data associated to padded array
    nghosts = PaddedArrays.npad(head_indices)  # number of ghost cells per boundary (compile-time constant)
    fill!(head_indices_data, EMPTY)  # this can be slow with too many cells?
    Np = length(xp)
    resize!(next_index, Np)
    Base.require_one_based_indexing(xp)
    Base.require_one_based_indexing(next_index)
    @inbounds Threads.@threads for n in 1:Np
        x⃗ = @inline get_coordinate(xp[n])  # usually get_coordinate === identity
        inds = map(determine_cell_index, Tuple(x⃗), rs_cell, Ls, size(cl))
        I = CartesianIndex(inds .+ nghosts)  # shift by number of ghost cells, since we access raw data associated to padded array
        head_old = Atomix.@atomicswap :monotonic head_indices_data[I] = n  # returns the old value
        # head_old = head_indices[I]
        # head_indices[I] = n       # the new element is the new head
        next_index[n] = head_old  # the old head now comes after the new element
    end
    pad_periodic!(cl.head_indices)  # fill ghost cells for periodicity (can be slow?)
    cl
end

set_elements!(cl::PeriodicCellList, xp::AbstractVector) = set_elements!(identity, cl, xp)

## ================================================================================ ##
## Iteration over a single cell

struct CellIterator{IndexType <: Integer}
    head_index :: IndexType          # index of first element in cell
    next_index :: Vector{IndexType}  # allows to get the rest of the elements in cell
end

Base.IteratorSize(::Type{<:CellIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:CellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:CellIterator{T}}) where {T} = T

function Base.iterate(it::CellIterator{IndexType}, n = it.head_index) where {IndexType}
    (; next_index,) = it
    n == IndexType(EMPTY) && return nothing  # no more elements in this cell
    @inbounds n, next_index[n]
end

## ================================================================================ ##
## Iteration over elements in cell lists

struct MultiCellIterator{
        IndexType <: Integer,
        N,
        HeadArray <: AbstractArray{IndexType, N},
        CellIndices <: CartesianIndices{N},
    }
    head_indices :: HeadArray          # array pointing to the index of the first element in each cell
    next_index   :: Vector{IndexType}  # allows to get the rest of the elements in cell
    cell_indices :: CellIndices        # indices of cells that will be visited
end

Base.IteratorSize(::Type{<:MultiCellIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:MultiCellIterator}) = Base.HasEltype()
Base.eltype(::Type{<:MultiCellIterator{IndexType}}) where {IndexType} = IndexType

@inline function Base.iterate(it::MultiCellIterator)
    icell = firstindex(it.cell_indices)
    I = @inbounds it.cell_indices[icell]
    n = @inbounds it.head_indices[I]
    iterate(it, (icell, n))
end

@inline function Base.iterate(it::MultiCellIterator{IndexType}, state) where {IndexType}
    (; head_indices, next_index, cell_indices,) = it
    icell, n = state
    while n == IndexType(EMPTY)
        if icell == lastindex(cell_indices)
            return nothing  # if this was the last cell, stop iterating
        end
        icell += 1  # jump to the next cell
        I = @inbounds cell_indices[icell]
        n = @inbounds head_indices[I]
    end
    @inbounds n, (icell, next_index[n])
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
    (; rs_cell, Ls,) = cl
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cell, Ls, size(cl))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    M = subdivisions(cl)
    cell_indices = CartesianIndices(
        map(i -> (i - M):(i + M), Tuple(I₀))
    )
    # iters = (CellIterator(cl.head_indices[I], cl.next_index) for I ∈ cell_indices)
    # it = Iterators.flatten(iters)  # this gives eltype(it) == Any (but that doesn't seem to affect performance...)
    it = MultiCellIterator(cl.head_indices, cl.next_index, cell_indices)  # might be slightly slower than Iterators.flatten (but with known eltype)
    it
end

end
