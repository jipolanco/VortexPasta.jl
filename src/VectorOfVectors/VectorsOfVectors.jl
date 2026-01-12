"""
    VectorsOfVectors

Module defining the `VectorOfVectors` type for conveniently working with values attached to filament nodes.
"""
module VectorsOfVectors

export VectorOfVectors

using Base: @propagate_inbounds

"""
    VectorOfVectors{T, V <: AbstractVector{T}} <: AbstractVector{V <: AbstractVector}
    VectorOfVectors(data::AbstractVector{<:AbstractVector})

Contains a list of vectors.

It behaves as much as possible as a basic vector of vectors. For instance:

- its `length` is the number of contained vectors;

- its element type (`eltype`) is the type of a contained vector;

- it only allows linear indexing (e.g. `u[i]` to get a single vector).
  To get an individual element (of type `T`), one should do `u[i][j]`.

Note that the individual vectors `u[i]` can have different lengths from each other.

There are also some differences:

- it overloads broadcasting mechanisms, so that doing `@. u = 2 * u + v` (where both
  variables are `VectorOfVectors`) works as expected and is efficient;

- copying a `VectorOfVectors` recursively copies its contained vectors, instead of just
  copying array references ("pointers").

# Examples

```jldoctest; setup = :(using Random; Random.seed!(42))
julia> data = [rand(n) for n ∈ 1:4]
4-element Vector{Vector{Float64}}:
 [0.6293451231426089]
 [0.4503389405961936, 0.47740714343281776]
 [0.7031298490032014, 0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278, 0.4570310908017041, 0.2993652953937611]

julia> us = VectorOfVectors(data)
4-element VectorOfVectors{Float64, Vector{Float64}}:
 [0.6293451231426089]
 [0.4503389405961936, 0.47740714343281776]
 [0.7031298490032014, 0.6733461456394962, 0.16589443479313404]
 [0.6134782250008441, 0.6683403279577278, 0.4570310908017041, 0.2993652953937611]

julia> us[2]
2-element Vector{Float64}:
 0.4503389405961936
 0.47740714343281776

julia> vs = @. us + 2 * us  # broadcasting
4-element VectorOfVectors{Float64, Vector{Float64}}:
 [1.888035369427827]
 [1.3510168217885807, 1.4322214302984533]
 [2.109389547009604, 2.0200384369184885, 0.4976833043794021]
 [1.8404346750025324, 2.0050209838731834, 1.3710932724051124, 0.8980958861812833]

```
"""
struct VectorOfVectors{T, V <: AbstractVector{T}} <: AbstractVector{V}
    u :: Vector{V}
end

Base.parent(x::VectorOfVectors) = x.u
Base.length(x::VectorOfVectors) = length(x.u)
Base.size(x::VectorOfVectors) = (length(x),)
Base.similar(x::VectorOfVectors, ::Type{T}) where {T} = map(v -> similar(v, T), x) :: VectorOfVectors
Base.similar(x::VectorOfVectors{T}) where {T} = similar(x, T)
Base.fill!(x::VectorOfVectors, value) = (foreach(y -> fill!(y, value), x); x)
Base.push!(x::VectorOfVectors, item) = push!(x.u, item)

Base.IndexStyle(::Type{<:VectorOfVectors}) = IndexLinear()
Base.IteratorSize(::Type{<:VectorOfVectors}) = Base.HasLength()

@propagate_inbounds Base.getindex(x::VectorOfVectors, i) = x.u[i]
@propagate_inbounds Base.setindex!(x::VectorOfVectors, v, i) = x.u[i] = v

Base.convert(::Type{VectorOfVectors}, data::AbstractVector) = VectorOfVectors(data)
Base.convert(::Type{VectorOfVectors}, x::VectorOfVectors) = x

# This is called when doing copy(us).
# The idea is to recursively copy each vector, kind of like `deepcopy`, instead of storing
# references to the same vectors.
function Base.copyto!(vs::VectorOfVectors, us::VectorOfVectors)
    axes(us) === axes(vs) || throw(DimensionMismatch("vectors have different sizes"))
    for (u, v) ∈ zip(us, vs)
        copyto!(v, u)
    end
    vs
end

# This is called when doing push!.
Base.resize!(vs::VectorOfVectors, n::Integer) = (resize!(vs.u, n); vs)
Base.pop!(vs::VectorOfVectors) = pop!(vs.u)
Base.popat!(vs::VectorOfVectors, i::Integer, args...) = popat!(vs.u, i, args...)
Base.empty!(vs::VectorOfVectors) = resize!(vs, 0)

# When generating a vector of vectors, make sure `map` returns a VectorOfVectors instead of
# a standard Vector.
function Base.map(f::F, vs::VectorOfVectors, args...) where {F}
    ys = map(f, vs.u, args...)
    _maybe_wrap(ys)
end

_maybe_wrap(ys::AbstractVector{<:AbstractVector}) = VectorOfVectors(ys)
_maybe_wrap(ys::AbstractVector{<:Any}) = ys

## Broadcasting

struct VectorOfVectorsStyle <: Broadcast.AbstractArrayStyle{1} end
VectorOfVectorsStyle(::Val{1}) = VectorOfVectorsStyle()

Base.BroadcastStyle(::Type{<:VectorOfVectors}) = VectorOfVectorsStyle()

function Base.similar(bc::Broadcast.Broadcasted{VectorOfVectorsStyle}, ::Type{ElType}) where {ElType}
    x = find_vov(bc)
    @assert ElType <: AbstractVector
    T = eltype(ElType)
    similar(x, T)
end

# This avoids obscure type instability when doing things like z = x .+ y.
# Note that the default definition is Broadcast.materialize(bc) = copy(instantiate(bc)),
# which is exactly the same thing we do below.
function Broadcast.materialize(bc::Broadcast.Broadcasted{VectorOfVectorsStyle})
    bc2 = Broadcast.instantiate(bc)
    copy(bc2)
end

# Find the first VectorOfVectors among the broadcasted arguments.
find_vov(bc::Broadcast.Broadcasted, args...) = find_vov(bc.args..., args...)
find_vov(u::VectorOfVectors, args...) = u
find_vov(::Any, args...) = find_vov(args...)
find_vov() = nothing  # this is never reached unless there's a bug somewhere...

# This custom implementation is there to make sure that broadcasting is done without
# unnecessary allocations.
function Base.copyto!(dest::VectorOfVectors, bc_in::Broadcast.Broadcasted{VectorOfVectorsStyle})
    for i ∈ eachindex(bc_in)
        bc = _extract_vectors(bc_in, i)
        Broadcast.materialize!(dest[i], bc)
    end
    dest
end

@inline function _extract_vectors(bc::Broadcast.Broadcasted, i)
    args = _extract_vectors(i, bc.args...)
    Base.broadcasted(bc.f, args...)
end

@inline _extract_vectors(i, arg, args...) =
    (_get_single_vector(arg, i), _extract_vectors(i, args...)...)
@inline _extract_vectors(i) = ()

@inline _get_single_vector(bc::Broadcast.Broadcasted, i) = _extract_vectors(bc, i)
@inline _get_single_vector(xs::AbstractVector, i) = @inbounds xs[i]  # includes VectorOfVectors
@inline _get_single_vector(x::Any, i) = x  # this can be the case of a scalar

end
