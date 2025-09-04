"""
    PaddedArrays

Module defining the [`PaddedArray`](@ref) type for dealing with arrays padded by ghost
cells.

Among other things, this type allows to efficiently work with periodic boundary conditions
in one and more spatial dimensions.
"""
module PaddedArrays

export PaddedArray, PaddedVector, pad_periodic!

"""
    PaddedArray{M, T, N} <: AbstractArray{T, N}

Pads a vector with `M` "ghost" entries on each side, along each direction.

Can be useful for dealing with periodic boundary conditions.
In that use case, one can call [`pad_periodic!`](@ref) once the non-ghost entries have been
filled to conveniently impose those kind of conditions.

See also [`PaddedVector`](@ref).

---

    PaddedArray{M}(data::AbstractArray)

Interpret input array as a padded array.

Note that the input array is not modified. Instead, its `M` first and `M` last
entries along each direction are considered as "ghost" entries.

In other words, the "logical" dimensions of the resulting `PaddedArray` are
`size(v) = size(data) .- 2M`.
Along a given direction of size `N`, indexing functions like `axes` return the range `1:N`
(or an equivalent).
However, the array can in reality be indexed (and modified) over the range `(1 - M):(N + M)`.

See [`PaddedVector`](@ref) for some one-dimensional examples.
"""
struct PaddedArray{M, T, N, A <: AbstractArray{T, N}} <: AbstractArray{T, N}
    data :: A
    function PaddedArray{M}(data::AbstractArray) where {M}
        A = typeof(data)
        T = eltype(A)
        N = ndims(A)
        Base.require_one_based_indexing(data)
        new{M, T, N, A}(data)
    end
end

npad(::Type{<:AbstractArray}) = 0  # regular arrays have no padding
npad(::Type{<:PaddedArray{M}}) where {M} = M
npad(v::AbstractArray) = npad(typeof(v))

# General case of N-D arrays: we can't use linear indexing due to the padding in all
# directions.
Base.IndexStyle(::Type{<:PaddedArray}) = IndexCartesian()

# Case of 1D arrays (vectors). Returns IndexLinear() if `A <: Vector`.
Base.IndexStyle(::Type{<:PaddedArray{M, T, 1, A}}) where {M, T, A} = IndexStyle(A)

Base.parent(v::PaddedArray) = v.data
Base.size(v::PaddedArray) = ntuple(i -> size(parent(v), i) - 2 * npad(v), Val(ndims(v)))

# Return a view to the parent array.
function Base.view(v::PaddedArray{M, T, N}, inds::Vararg{Any, N}) where {M, T, N}
    inds_parent = map(inds, axes(v)) do is, js
        to_parent_indices(M, is, js)
    end
    view(parent(v), inds_parent...)
end

to_parent_indices(M, i::Integer, js) = M + i
to_parent_indices(M, is::AbstractVector, js) = M .+ is  # this includes ranges
to_parent_indices(M, ::Colon, js) = M .+ js  # select all values excluding ghost cells

function Base.copyto!(w::PaddedArray{M}, v::PaddedArray{M}) where {M}
    # If v and w point to the same array, there's no need to copy data.
    # This also avoids error when data arrays are UnsafeArrays returned by Bumper.
    w.data === v.data && return w
    length(w.data) == length(v.data) || throw(DimensionMismatch("arrays have different sizes"))
    copyto!(w.data, v.data)
    w
end

# Equality and isapprox must also include ghost cells!
Base.:(==)(w::PaddedArray{M}, v::PaddedArray{M}) where {M} = w.data == v.data

Base.isapprox(w::PaddedArray{M}, v::PaddedArray{M}; kwargs...) where {M} =
    isapprox(w.data, v.data; kwargs...)

Base.fill!(v::PaddedArray, x) = fill!(parent(v), x)

function Base.similar(v::PaddedArray{M, T, N}, ::Type{S}, dims::Dims{N}) where {S, M, T, N}
    data = similar(v.data, S, dims .+ 2M)
    fill!(data, zero(S))
    PaddedArray{M}(data)
end

Base.@propagate_inbounds function Base.getindex(v::PaddedArray{M, T, N}, I::Vararg{Int, N}) where {M, T, N}
    J = ntuple(n -> I[n] + M, Val(N))
    parent(v)[J...]
end

Base.@propagate_inbounds function Base.setindex!(
        v::PaddedArray{M, T, N},
        val,
        I::Vararg{Int, N},
    ) where {M, T, N}
    J = I .+ M
    parent(v)[J...] = val
end

Base.checkbounds(::Type{Bool}, v::PaddedArray, I...) = _checkbounds(v, I...)

_checkbounds(v::PaddedArray, I::CartesianIndex) = _checkbounds(v, Tuple(I)...)

function _checkbounds(v::PaddedArray, Is...)
    M = npad(v)
    N = ndims(v)
    P = length(Is)
    @assert P ≥ N
    Is_head = ntuple(n -> Is[n], Val(N))
    Is_tail = ntuple(n -> Is[N + n], Val(P - N))  # possible additional indices which should be all 1
    for (inds, i) ∈ zip(axes(v), Is_head)
        _checkbounds(Val(M), inds, i) || return false
    end
    all(isone, Is_tail)
end

_checkbounds(::Val{M}, inds::AbstractUnitRange, i::Integer) where {M} =
    first(inds) - M ≤ i ≤ last(inds) + M

_checkbounds(::Val{M}, inds::AbstractUnitRange, I::AbstractUnitRange) where {M} =
    first(inds) - M ≤ first(I) && last(I) ≤ last(inds) + M

## ================================================================================ ##
## Specialisations for the 1D case (PaddedVector).
## ================================================================================ ##

"""
    PaddedVector{M, T} <: AbstractVector{T}

Alias for `PaddedArray{M, T, 1}` which can be used to work with one-dimensional data.

---

    PaddedVector{M}(data::AbstractVector)

Interpret input vector as a padded vector.

See [`PaddedArray`](@ref) for details.

# Examples

```jldoctest
julia> v = PaddedVector{2}(collect(1:10))
6-element PaddedVector{2, Int64, Vector{Int64}}:
 3
 4
 5
 6
 7
 8

julia> eachindex(v)
Base.OneTo(6)

julia> v[begin]
3

julia> v[begin - 2]
1

julia> v[end]
8

julia> v[end + 2]
10

julia> v[end - 1] = 42; println(v)
[3, 4, 5, 6, 42, 8]
```
"""
const PaddedVector{M, T, V} = PaddedArray{M, T, 1, V}

PaddedVector{M}(data::AbstractVector) where {M} = PaddedArray{M}(data)

function Base.resize!(v::PaddedVector, n::Integer)
    n_prev = length(v)
    if n_prev != n
        resize!(parent(v), n + 2 * npad(v))
        fill!(parent(v), zero(eltype(v)))  # reset all values just in case
    end
    v
end

function Base.sizehint!(v::PaddedVector, n::Integer)
    sizehint!(parent(v), n + 2 * npad(v))
    v
end

function Base.insert!(v::PaddedVector, i::Integer, x)
    insert!(parent(v), npad(v) + i, x)
    v
end

function Base.push!(v::PaddedVector, x)
    resize!(v, length(v) + 1)
    v[end] = x
    v
end

Base.popat!(v::PaddedVector, i::Integer) = popat!(parent(v), i + npad(v))

# This is used by HDF5.jl when reading data directly onto a PaddedVector using HDF5.API.h5d_read()
Base.unsafe_convert(::Type{Ptr{T}}, v::PaddedVector{M, T}) where {M, T} = pointer(v)
Base.pointer(v::PaddedVector) = pointer(parent(v), npad(v) + 1)  # points to the first non-ghost entry

## ================================================================================ ##
## Periodic padding.
## ================================================================================ ##

struct FromCentre end
struct FromRight end
struct FromLeft end

"""
    pad_periodic!(v::PaddedArray{M, T, N}, [L = zero(T)])

Fill ghost cells in a periodic manner.

In the simplest case of a 1D `PaddedArray` (`N = 1`), this function will copy:

- `v[begin:(begin + M - 1)]` → `v[(end + 1):(end + M)]`, and

- `v[(end - M + 1):end]` → `v[(begin - M):(begin - 1)]`.

Something equivalent (but more complicated) is done in multiple dimensions.

If `L ≠ 0`, it is interpreted as an unfolding period (or as an end-to-end distance after a
single period), such that `v[N + i] - v[i] = L`, where `N = length(v)`.
"""
function pad_periodic! end

pad_periodic!(v::PaddedArray, args...) = pad_periodic!(FromCentre(), v, args...)

function pad_periodic!(::FromCentre, v::PaddedVector{M, T}, L::T = zero(T)) where {M, T}
    @inbounds for i ∈ 1:M
        v[begin - i] = v[end + 1 - i] - L
        v[end + i] = v[begin - 1 + i] + L
    end
    v
end

# This variant gives priority to padded values on the right of the "central" array.
# This can be convenient for certain algorithms (e.g. when inserting spline knots).
function pad_periodic!(::FromRight, v::PaddedVector{M, T}, L::T = zero(T)) where {M, T}
    # Copy right ghost values onto "central" array, then perform the usual padding "from centre".
    @inbounds for i ∈ 1:M
        v[begin - 1 + i] = v[end + i] - L
    end
    pad_periodic!(FromCentre(), v, L)
end

# Similar to above, but preserving the left padding.
function pad_periodic!(::FromLeft, v::PaddedVector{M, T}, L::T = zero(T)) where {M, T}
    # Copy left ghost values onto "central" array, then perform the usual padding "from centre".
    @inbounds for i ∈ 1:M
        v[end + 1 - i] = v[begin - i] + L
    end
    pad_periodic!(FromCentre(), v, L)
end

# Generalisation to N-dimensional PaddedArray.
@inline function pad_periodic!(dir::FromCentre, v::PaddedArray)
    # Copy ghost cells over each dimension separately.
    N = ndims(v)
    _pad_periodic_dimension!(dir, Val(N), v)
end

# Copy values over dimension d.
@inline function _pad_periodic_dimension!(dir::FromCentre, ::Val{d}, v::PaddedArray) where {d}
    data = parent(v)
    N = ndims(data)
    inds_before = CartesianIndices(ntuple(i -> axes(data, i), Val(d - 1)))
    inds_after = CartesianIndices(ntuple(i -> axes(data, d + i), Val(N - d)))
    ibegin = firstindex(data, d)
    iend = lastindex(data, d)
    M = npad(v)
    @inbounds for J ∈ inds_after, I ∈ inds_before
        for δ ∈ 1:M  # copy each ghost cell layer
            data[I, iend - M + δ, J] = data[I, ibegin - 1 + M + δ, J]
            data[I, ibegin - 1 + δ, J] = data[I, iend - 2M + δ, J]
        end
    end
    _pad_periodic_dimension!(dir, Val(d - 1), v)
end

# Stop when we have copied over all dimensions.
@inline _pad_periodic_dimension!(dir::FromCentre, ::Val{0}, v::PaddedArray) = v

end  # module
