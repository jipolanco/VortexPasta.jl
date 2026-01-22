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
public set_elements!, foreach_pair, foreach_source

using ..PaddedArrays: PaddedArrays, PaddedArray, pad_periodic!
using Atomix: Atomix
using KernelAbstractions: KernelAbstractions as KA, CPU, GPU
using Static: StaticInt, static, dynamic
using StaticArrays: MVector

# This is used either to mean that a cell has no elements, or that an element is the last
# element within a given cell. This assumes one-based indexing, so that 0 is an invalid index!
const EMPTY = 0

## ================================================================================ ##
## Cell list type definition

"""
    PeriodicCellList(
        [backend = CPU()],
        rs_cut::NTuple{N, Real}, periods::NTuple{N, Real},
        [nsubdiv = Val(1)];
    ) -> PeriodicCellList{N}

Construct a cell list for dealing with pair interactions in `N` dimensions.

The cutoff distances `rs_cut` (which can be different in each direction) don't need to exactly
divide the domain period `L` into equal pieces. In that case the cutoff distances will be
adjusted (slightly increased) to fit the domain period.

Note that this cell-lists implementation can identify pairs which are slightly beyond the
given cut-off distance. For performance reasons, such pairs are still returned, and it's up
to the user to see whether such pairs should still be computed or not depending on whether a
strict cut-off distance is wanted.

Optionally, one can choose to subdivide each cell (of size `≈ rcut`) onto `nsubdiv`
subcells. This can significantly improve performance, since it allows to discard some
spurious pair interactions (i.e. beyond the chosen cutoff radius) as described
[here](https://en.wikipedia.org/wiki/Cell_lists#Improvements). In practice, a value of
`2` or `3` can significantly improve performance compared to no subdivision (`1`).
For convenience, it can be passed as `nsubdiv = Val(M)` or as `nsubdiv = static(M)`.

Infinite non-periodic domains (in the sense of `period = Infinity()`) are not supported.

## GPU usage

To run on GPUs, pass a KernelAbstractions backend such as `CUDABackend` or `ROCBackend` as
the first argument.

If multiple GPUs are available, make sure to activate the device where computations should
be performed _before_ constructing a `PeriodicCellList`. For example:

```julia
using KernelAbstractions: KernelAbstractions as KA
backend = ROCBackend()
device_id = 2  # assuming there are two or more available GPUs
KA.device!(backend, device_id)
cl = PeriodicCellList(backend, ...)
```
"""
struct PeriodicCellList{
        N,  # spatial dimension
        M,  # number of cell subdivisions (≥ 1)
        KABackend <: KA.Backend,
        HeadArray <: PaddedArray{M, Int, N},  # array of cells, each containing a list of elements
        IndexVector <: AbstractVector{Int},
        CellDims <: NTuple{N, Real},
        Periods <: NTuple{N, Real},
    }
    backend       :: KABackend    # CPU, CUDABackend, ROCBackend, ...
    device        :: Int          # device id (in 1:ndevices) where arrays are stored
    head_indices  :: HeadArray    # array pointing to the index of the first element in each cell (or EMPTY if no elements in that cell)
    next_index    :: IndexVector  # points from one element to the next element (or EMPTY if this is the last element) [length Np]
    rs_cell       :: CellDims     # dimensions of a cell (can be different in each direction)
    nsubdiv       :: StaticInt{M}
    Ls            :: Periods
end

function Base.show(io::IO, cl::PeriodicCellList{N}) where {N}
    (; backend, rs_cell, Ls, nsubdiv) = cl
    print(io, "PeriodicCellList{$N} with:")
    print(io, "\n - backend:         ", backend)
    print(io, "\n - number of cells: ", size(cl))
    print(io, "\n - period:          ", Ls)
    print(io, "\n - cell size:       ", rs_cell)
    print(io, "\n - effective r_cut: ", rs_cell .* dynamic(nsubdiv))
    print(io, "\n - number of subdivisions: ", dynamic(nsubdiv))
    nothing
end

Base.size(cl::PeriodicCellList) = size(cl.head_indices)
subdivisions(::PeriodicCellList{N, M}) where {N, M} = M

@inline PeriodicCellList(rs_cut::NTuple{N, Real}, args...; kws...) where {N} = PeriodicCellList(CPU(), rs_cut, args...; kws...)
@inline PeriodicCellList(backend, rs, Ls, nsubdiv::Val{M}; kws...) where {M} = PeriodicCellList(backend, rs, Ls, static(M); kws...)

# Returns maximum possible cut-off distance r_cut for M subdivisions and a domain of period L.
max_cutoff_distance(M::Integer, L::AbstractFloat) = oftype(L, M / (2M + 1)) * L

function PeriodicCellList(
        backend::KA.Backend,
        rs_cut::NTuple{N, Real},
        Ls::NTuple{N, Real},
        nsubdiv::StaticInt = static(1),
    ) where {N}
    any(isinf, Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by PeriodicCellList"
    ))

    device = KA.device(backend)  # currently active device

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
    head_raw = KA.allocate(backend, IndexType, head_dims)
    fill!(head_raw, EMPTY)

    head_indices = PaddedArray{M}(head_raw)
    @assert size(head_indices) == ncells

    next_index = KA.allocate(backend, IndexType, 0)

    PeriodicCellList(backend, device, head_indices, next_index, rs_cell, nsubdiv, Ls)
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
    determine_cell_index_folded(x, rcut, N)
end

# Like determine_cell_index, but assumes that `x` is already in [0, L).
@inline function determine_cell_index_folded(x, rcut, N)
    # Here unsafe_trunc(Int, ⋅) is used instead of floor(Int, ⋅) because it should be faster.
    # For non-negative values, both should give the same result.
    # The unsafe_trunc function generally uses a single intrinsic CPU instruction and never
    # throws errors. It can silently give a wrong result if the values are not representable
    # by an Int, but that will never be the case in practice here (since 0 ≤ x/rcut < L/rcut
    # and L/rcut is very small compared to typemax(Int) = 2^63 - 1).
    clamp(1 + unsafe_trunc(Int, x / rcut), 1, N)  # make sure the index is in 1:N
end

@inline function get_neighbouring_cell_indices(cl::PeriodicCellList, x⃗, folded::Val = Val(false))
    (; rs_cell, Ls) = cl
    if folded === Val(true)
        inds_central = map(determine_cell_index_folded, Tuple(x⃗), rs_cell, size(cl))
    else
        inds_central = map(determine_cell_index, Tuple(x⃗), rs_cell, Ls, size(cl))
    end
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    M = subdivisions(cl)
    CartesianIndices(map(i -> (i - M):(i + M), Tuple(I₀)))
end

"""
    CellLists.set_elements!([get_coordinate::Function], cl::PeriodicCellList, xp::AbstractVector; folded = Val(false))

Set all elements of the cell list.

In general `xp` is a vector of spatial locations.

Optionally, a `get_coordinate` function can be passed, which will be applied to each `xp[j]` in order to obtain a coordinate.
By default this is `identity`.

Moreover, one may pass `folded = Val(true)` in case all points in `xp` have been folded into
the main periodic cell (i.e. `x⃗ ∈ [0, L]ᵈ`), which can speed-up computations. **This is a
dangerous option**: if points haven't been really folded, this lead to wrong results or crashes.

This function resets the cell list, removing all previously existent points.
"""
function set_elements!(
        get_coordinate::F, cl::PeriodicCellList, xp::AbstractVector;
        folded::Val{Folded} = Val(false),
    ) where {F, Folded}
    Folded::Bool
    (; backend, device) = cl
    # typeof(KA.get_backend(xp)) === typeof(backend) ||
    #     throw(ArgumentError(lazy"coordinate vector `xp` should be on the computing device ($backend)"))
    device_current = KA.device(backend)
    device == device_current ||
        error(lazy"the currently selected device (id = $device_current) is not the same that was used to create the PeriodicCellList (id = $device)")
    _set_elements!(backend, get_coordinate, cl, xp, folded)
end

set_elements!(cl::PeriodicCellList, xp::AbstractVector; kws...) = set_elements!(identity, cl, xp; kws...)

function _set_elements!(::CPU, get_coordinate::F, cl::PeriodicCellList, xp::AbstractVector, folded::Val) where {F}
    (; next_index, head_indices, Ls, rs_cell,) = cl
    head_indices_data = parent(head_indices)   # full data associated to padded array
    nghosts = PaddedArrays.npad(head_indices)  # number of ghost cells per boundary (compile-time constant)
    Threads.@threads for i in eachindex(head_indices_data)
        @inbounds head_indices_data[i] = EMPTY
    end
    Np = length(xp)
    resize!(next_index, Np)
    Base.require_one_based_indexing(xp)
    Base.require_one_based_indexing(next_index)
    @inbounds Threads.@threads for n in 1:Np
        x⃗ = @inline get_coordinate(xp[n])  # usually get_coordinate === identity
        if folded === Val(true)
            inds = map(determine_cell_index_folded, Tuple(x⃗), rs_cell, size(cl))
        else
            inds = map(determine_cell_index, Tuple(x⃗), rs_cell, Ls, size(cl))
        end
        I = CartesianIndex(inds .+ nghosts)  # shift by number of ghost cells, since we access raw data associated to padded array
        head_old = Atomix.@atomicswap :monotonic head_indices_data[I] = n  # returns the old value
        # head_old = head_indices[I]
        # head_indices[I] = n       # the new element is the new head
        next_index[n] = head_old  # the old head now comes after the new element
    end
    pad_periodic!(cl.head_indices)  # fill ghost cells for periodicity
    cl
end

"""
    CellLists.foreach_pair(
        f::Function, cl::PeriodicCellList, xp_dest::AbstractVector;
        batch_size = nothing,
        folded = Val(false),
    )

Iterate over all point pairs within the chosen cut-off distance.

A pair is given by a **destination point** `xp_dest[i]` and a **source point** `xp[j]`,
where `xp` is the vector of coordinates that was previously passed to [`CellLists.set_elements!`](@ref).
Note that the vector of destination points `xp_dest` can be different than that of source points.

The first argument should be a user-defined function `f(x⃗, i::Integer, j::Integer)` taking a
destination point `x⃗ = xp_dest[i]` and a pair of destination/source indices `(i, j)`.
Typically, one wants to modify one or more destination vectors at index `vp[i]`, which can
be done within this function.

On CPUs this function is parallelised using threads. Parallelisation is done at the level of
all destination points `xp_dest`, ensuring that only a single thread can write to a location
`vp[i]` (in other words, atomics are not needed). This function also works on GPUs.

One may set `folded = Val(true)` if points in `xp_dest` have been previously folded.
See [`set_elements!`](@ref) for more details.

This function should be called after [`CellLists.set_elements!`](@ref).

See also [`CellLists.foreach_source`](@ref).

# Iterating over batches

Similarly to [`foreach_source`](@ref), one can iterate over batches of (up to) `W` source
points at a time (for each single destination point `x⃗ = xp_dest[i]`). For this one should
pass `batch_size = Val(W)` as a keyword argument. Note that in this case, the user-defined
function `f` is expected to have the signature `f(x⃗, i::Integer, js::NTuple{W}, m::Integer)`.
See [`foreach_source`](@ref) for details on the last two arguments.

# Example usage

```jldoctest
julia> using StaticArrays: SVector

julia> using LinearAlgebra: ⋅  # dot product

julia> r_cut = 0.4; rs_cut = (r_cut, r_cut, r_cut);

julia> Ls = (2π, 2π, 2π);

julia> nsubdiv = Val(2);

julia> cl = PeriodicCellList(CPU(), rs_cut, Ls, nsubdiv)
PeriodicCellList{3} with:
 - backend:         CPU(false)
 - number of cells: (31, 31, 31)
 - period:          (6.283185307179586, 6.283185307179586, 6.283185307179586)
 - cell size:       (0.2026833970057931, 0.2026833970057931, 0.2026833970057931)
 - effective r_cut: (0.4053667940115862, 0.4053667940115862, 0.4053667940115862)
 - number of subdivisions: 2

julia> xp = [rand(SVector{3, Float64}) .* Ls for _ in 1:1000];  # create source points in [0, L]³

julia> CellLists.set_elements!(cl, xp);  # process source points

julia> xp_dest = [rand(SVector{3, Float64}) .* Ls for _ in 1:200];  # create destination points in [0, L]³

julia> vp_dest = zeros(length(xp_dest));  # vector where "interactions" will be written to

julia> CellLists.foreach_pair(cl, xp_dest) do x⃗, i, j
           # Compute some "influence" of xp[j] on x⃗ == xp_dest[i]
           # Note that it is safe to write to vp_dest[i] without atomics, even when running in parallel.
           @inbounds vp_dest[i] += x⃗ ⋅ xp[j]
       end
```
"""
function foreach_pair(
        f::F, cl::PeriodicCellList, xp_dest::AbstractVector;
        batch_size::Union{Nothing, Val} = nothing,
        folded::Val{Folded} = Val(false),
    ) where {F <: Function, Folded}
    Folded::Bool
    (; backend) = cl
    _foreach_pair(backend, f, cl, xp_dest, batch_size, folded)
    nothing
end

# This is similar to Base.Fix / Fix1 / Fix2.
struct Fix12{F <: Function, A, B} <: Function
    f::F
    a::A
    b::B
end

@inline (g::Fix12)(args::Vararg{Any, N}) where {N} = @inline g.f(g.a, g.b, args...)

# On the CPU, this simply calls foreach_source in parallel for each destination point.
function _foreach_pair(backend::CPU, f::F, cl, xp_dest, batch_size, folded::Val) where {F}
    Threads.@threads :dynamic for i in eachindex(xp_dest)
        @inbounds x⃗ = xp_dest[i]
        g = Fix12(f, x⃗, i)
        _foreach_source(backend, g, cl, x⃗, batch_size, folded)
    end
    nothing
end

"""
    CellLists.foreach_source(f::Function, cl::PeriodicCellList, x⃗_dest; batch_size = nothing, folded = Val(false))

Iterate over source points which are close to a given destination point `x⃗_dest`.

This is similar to [`CellLists.foreach_pair`](@ref) but only works on a single destination
point. On the CPU, `foreach_source` is not parallelised. One may want to parallelise at the
level of multiple destination points to match the behaviour of `foreach_pair`.

The first argument should be a user-defined function `f(j::Integer)` taking the index `j` of
a source point.

Performance-wise, this function seems to give similar performance to `foreach_pair` on CPUs.
On GPUs one may prefer to use `foreach_pair` which is more adapted for GPU computing.

One may set `folded = Val(true)` if `x⃗_dest` has been previously folded.
See [`set_elements!`](@ref) for more details.

# Iterating over batches

For performance reasons (e.g. to enable SIMD), one may want to iterate over batches of
source points in `cl`. For this, one should pass the keyword argument `batch_size = Val(W)`
where `W` is the wanted batch size. 

When batching is enabled, the user-defined function `f` is expected to have the signature
`f(js::NTuple{W}, m::Integer)`, where `js` are the indices of `W` nearby source points, and
`m ∈ 1:W` is the number of "valid" points. In other words, indices `js[(m + 1):W]` should be
discarded (or masked out) in the `f` function. Note that, even though these indices should
be ignored to get correct results, these are guaranteed to be valid indices in `1:Np` (where
`Np` is the number of source points), and thus can be safely used to access source
coordinates (which is convenient in the context of SIMD).
"""
@inline function foreach_source(
        f::F, cl::PeriodicCellList, x⃑_dest;
        batch_size::Union{Nothing, Val} = nothing,
        folded::Val{Folded} = Val(false),
    ) where {F <: Function, Folded}
    Folded::Bool
    (; backend) = cl
    _foreach_source(backend, f, cl, x⃑_dest, batch_size, folded)
    nothing
end

@inline function _foreach_source(::CPU, f::F, cl, x⃗, batch_size::Nothing, folded::Val) where {F}
    (; head_indices, next_index) = cl
    cell_indices = get_neighbouring_cell_indices(cl, x⃗, folded)
    # Iterate over "active" cells around the destination point.
    # TODO: exclude "corner" cells which are for sure beyond the cut-off distance?
    @inbounds for I in cell_indices
        j = head_indices[I]
        while j != EMPTY
            @inline f(j)
            j = next_index[j]
        end
    end
    nothing
end

@inline function _foreach_source(::CPU, f::F, cl, x⃗, ::Val{batch_size}, folded::Val) where {F, batch_size}
    (; head_indices, next_index) = cl
    IndexType = eltype(head_indices)
    cell_indices = get_neighbouring_cell_indices(cl, x⃗, folded)
    inds = MVector{batch_size, IndexType}(undef)
    m = 0
    @inbounds for I in cell_indices
        j = head_indices[I]
        while j != EMPTY
            inds[m += 1] = j
            if m == batch_size
                @inline f(Tuple(inds), batch_size)
                m = 0
            end
            j = next_index[j]
        end
    end
    if m > 0
        for l in (m + 1):batch_size
            # copy latest index, just to make sure that all returned indices are valid
            @inbounds inds[l] = inds[m]
        end
        @inline f(Tuple(inds), m)
    end
    nothing
end

## ================================================================================ ##
## Iterator interface: iteration over elements in cell lists

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
    cell_indices = get_neighbouring_cell_indices(cl, x⃗)
    # iters = (CellIterator(cl.head_indices[I], cl.next_index) for I ∈ cell_indices)
    # it = Iterators.flatten(iters)  # this gives eltype(it) == Any (but that doesn't seem to affect performance...)
    it = MultiCellIterator(cl.head_indices, cl.next_index, cell_indices)  # might be slightly slower than Iterators.flatten (but with known eltype)
    it
end

include("gpu.jl")

end
