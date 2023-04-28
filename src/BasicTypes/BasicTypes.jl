"""
    BasicTypes

Defines and exports basic types used in all computations.
"""
module BasicTypes

export
    Vec3,
    Derivative,
    Zero,
    Infinity

using StaticArrays: SVector

"""
    Vec3{T}

Three-element static vector, alias for `SVector{3, T}`.

Used to describe vectors and coordinates in 3D space.
"""
const Vec3{T} = SVector{3, T}
    
"""
    Derivative{N}

Represents the ``N``-th order derivative operator.

Used in particular to interpolate derivatives along filaments.
"""
struct Derivative{N} end
Derivative(N::Int) = Derivative{N}()
Base.broadcastable(d::Derivative) = Ref(d)  # disable broadcasting on Derivative objects

include("constants.jl")

end
