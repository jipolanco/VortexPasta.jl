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

# Used internally to evaluate filament coordinates or derivatives on a given
# discretisation node. This is used when calling f[i, Derivative(n)].
struct AtNode
    i :: Int
end
