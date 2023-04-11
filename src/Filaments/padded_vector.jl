export PaddedVector

"""
    PaddedVector{M, T} <: AbstractVector{T}

Pads a vector with `M` "ghost" entries on each side.

Can be useful for dealing with periodic boundary conditions (in our case for
closed filaments).

---

    PaddedVector{M}(data::AbstractVector)

Interpret input array as a padded vector.

Note that the input array is not modified. Instead, its `M` first and `M` last
entries are considered as "ghost" entries.

In other words, the "logical" length of the resulting `PaddedVector` is
`N = length(data) - 2M`.
Functions like `eachindex` or `axes` return the range `1:N` (or an equivalent).
However, the padded vector can in reality be indexed (and modified) over the range
`(1 - M):(N + M)`.

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
struct PaddedVector{M, T, V <: AbstractVector{T}} <: AbstractVector{T}
    data :: V
    function PaddedVector{M}(data::AbstractVector) where {M}
        T = eltype(data)
        V = typeof(data)
        Base.require_one_based_indexing(data)
        new{M, T, V}(data)
    end
end

npad(::Type{<:PaddedVector{M}}) where {M} = M
npad(v::PaddedVector) = npad(typeof(v))

# Not sure if this makes a difference for 1D arrays...
Base.IndexStyle(::Type{<:PaddedVector{M, T, V}}) where {M, T, V} = IndexStyle(V)

Base.parent(v::PaddedVector) = v.data
Base.size(v::PaddedVector) = (length(parent(v)) - 2 * npad(v),)
Base.size(v::PaddedVector{0}) = size(parent(v))

Base.checkbounds(::Type{Bool}, v::PaddedVector, I...) = _checkbounds(v, I...)

# TODO define copyto! for different paddings M, M′?
function Base.copyto!(w::PaddedVector{M}, v::PaddedVector{M}) where {M}
    copyto!(w.data, v.data)
    w
end

@inline function _checkbounds(v::PaddedVector, i::Integer)
    M = npad(v)
    firstindex(v) - M ≤ i ≤ lastindex(v) + M
end

@inline function _checkbounds(v::PaddedVector, I::AbstractUnitRange)
    M = npad(v)
    firstindex(v) - M ≤ first(I) && last(I) ≤ lastindex(v) + M
end

# For some reason, when printing `v`, Julia indexes it as a 2D array (with j = 1).
# This is needed to avoid failure.
@inline function _checkbounds(v::PaddedVector, i::Integer, Is::Vararg{Integer})
    all(isone, Is) && _checkbounds(v, i)
end

function Base.similar(v::PaddedVector, ::Type{S}, dims::Dims{1}) where {S}
    M = npad(v)
    PaddedVector{M}(similar(v.data, S, dims .+ 2M))
end

Base.@propagate_inbounds Base.getindex(v::PaddedVector, i::Int) =
    parent(v)[i + npad(v)]

Base.@propagate_inbounds Base.setindex!(v::PaddedVector, val, i::Int) =
    parent(v)[i + npad(v)] = val

# Apply periodic padding.
# If L ≠ 0, it is interpreted as an unfolding period, such that v[N + 1 + i] - v[i] = L (where N = length(v)).
function pad_periodic!(v::PaddedVector{M, T}, L::T = zero(T)) where {M, T}
    if length(v) ≥ M
        @inbounds for i ∈ 1:M
            v[begin - i] = v[end + 1 - i] - L
            v[end + i] = v[begin - 1 + i] + L
        end
    else
        # TODO apply "partial" (or multiple?) periodic padding based on the
        # number of elements
    end
    v
end
