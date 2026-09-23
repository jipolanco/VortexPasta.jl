using StaticArrays: SVector

"""
    Vec3{T}

Three-element static vector, alias for `SVector{3, T}`.

Used to describe vectors and coordinates in 3D space.
"""
const Vec3{T} = SVector{3, T}

# This is copied from BSplineKit.jl.
abstract type AbstractDifferentialOp end
    
"""
    Derivative{N}

Represents the ``N``-th order derivative operator.

Used in particular to interpolate derivatives along filaments.
"""
struct Derivative{N} <: AbstractDifferentialOp end
Derivative(N::Int) = Derivative{N}()
Base.broadcastable(d::Derivative) = Ref(d)  # disable broadcasting on Derivative objects
output_type(::Derivative, ::Type{T}) where {T} = T

# This is for internal use only (splines) -- adapted from BSplineKit.jl
struct ScaledDifferentialOp{Op <: AbstractDifferentialOp, T <: Number} <: AbstractDifferentialOp
    D::Op
    α::T
end

output_type(S::ScaledDifferentialOp, ::Type{T}) where {T} = promote_type(output_type(S.D, T), typeof(S.α))
Base.show(io::IO, S::ScaledDifferentialOp) = print(io, "(", S.α, ") * ", S.D)
Base.:-(S::ScaledDifferentialOp) = ScaledDifferentialOp(S.D, -S.α)
Base.:*(α::Number, D::AbstractDifferentialOp) = ScaledDifferentialOp(D, α)
Base.:*(D::AbstractDifferentialOp, α) = α * D
Base.:-(D::AbstractDifferentialOp) = -1 * D

# This is for internal use only (splines) -- adapted from BSplineKit.jl
struct DifferentialOpSum{
        OpA <: AbstractDifferentialOp,
        OpB <: AbstractDifferentialOp,
    } <: AbstractDifferentialOp
    a::OpA
    b::OpB
end

output_type(op::DifferentialOpSum, ::Type{T}) where {T} = promote_type(output_type(op.a, T), output_type(op.b, T))
Base.show(io::IO, op::DifferentialOpSum) = print(io, op.a, " + ", op.b)
Base.:+(a::AbstractDifferentialOp, b::AbstractDifferentialOp) = DifferentialOpSum(a, b)
Base.:-(a::AbstractDifferentialOp, b::AbstractDifferentialOp) = DifferentialOpSum(a, -b)

# Used internally to evaluate filament coordinates or derivatives on a given
# discretisation node. This is used when calling f[i, Derivative(n)].
struct AtNode
    i :: Int
end
