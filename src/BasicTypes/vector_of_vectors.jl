using Base: @propagate_inbounds

"""
    VectorOfVectors{T, V <: AbstractVector{T}} <: AbstractVector{V <: AbstractVector}

Contains a list of vectors.

It behaves as much as possible as a vector of vectors, with some convenience additions e.g.
for broadcasting over vectors of vectors.

This is somewhat similar to the `VectorOfArray` type from the `RecursiveArrayTools.jl`
package, but with some important differences:

- a `VectorOfVectors` only allows linear indexing (e.g. `u[i]` to get a single vector).
  To get an individual element (of type `T`), one should do `u[i][j]`.

- a `VectorOfVectors` has `IndexStyle(u) = IndexLinear()`, which means that
  `eachindex(u, v)` returns a linear range.

- the element type is an `AbstractVector{T}` instead of `T`.

- defines functions such as `popat!`.
"""
struct VectorOfVectors{T, V <: AbstractVector{T}} <: AbstractVector{V}
    u :: Vector{V}
end

Base.length(x::VectorOfVectors) = length(x.u)
Base.size(x::VectorOfVectors) = (length(x),)
Base.similar(x::VectorOfVectors, ::Type{T}) where {T} =
    VectorOfVectors(map(v -> similar(v, T), x.u))
Base.similar(x::VectorOfVectors{T}) where {T} = similar(x, T)

Base.IndexStyle(::Type{<:VectorOfVectors}) = IndexLinear()
Base.IteratorSize(::Type{<:VectorOfVectors}) = Base.HasLength()

@propagate_inbounds Base.getindex(x::VectorOfVectors, i) = x.u[i]
@propagate_inbounds Base.setindex!(x::VectorOfVectors, v, i) = x.u[i] = v

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
Base.resize!(vs::VectorOfVectors, n::Integer) = resize!(vs.u, n)

Base.pop!(vs::VectorOfVectors) = pop!(vs.u)
Base.popat!(vs::VectorOfVectors, i::Integer, args...) = popat!(vs.u, i, args...)

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

# Find the first VectorOfVectors among the broadcasted arguments.
find_vov(bc::Broadcast.Broadcasted) = find_vov(bc.args...)
find_vov(u::VectorOfVectors, args...) = u
find_vov(::Any, args...) = find_vov(args...)

function Base.copyto!(dest::VectorOfVectors, bc_in::Broadcast.Broadcasted{VectorOfVectorsStyle})
    bc = Broadcast.flatten(bc_in)
    for i ∈ eachindex(bc)
        # This is to make sure that broadcasting is done without unnecessary allocations.
        args = map(w -> _get_single_vector(w, i), bc.args)
        broadcast!(bc.f, dest[i], args...)
    end
    dest
end

_get_single_vector(x::VectorOfVectors, i) = @inbounds x[i]
_get_single_vector(x::Any, i) = x  # this can be the case of a scalar
