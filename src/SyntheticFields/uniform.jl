"""
    UniformVectorField{T, N} <: SyntheticVectorField{T, N}
    UniformVectorField(u⃗)

Represents a uniform (constant) vector field with no spatial fluctuations.

This can be used for instance to represent a uniform counterflow.

The input `u⃗` can be an `SVector{N, T}` or a tuple.
"""
struct UniformVectorField{T, N} <: SyntheticVectorField{T, N}
    u⃗::SVector{N, T}
end

function Base.show(io::IO, field::UniformVectorField{T, N}) where {T, N}
    print(io, "UniformVectorField{$T, $N} with value ", field.u⃗)
end

UniformVectorField(u⃗::Tuple) = UniformVectorField(SVector(u⃗))

(field::UniformVectorField)(x⃗) = field.u⃗
