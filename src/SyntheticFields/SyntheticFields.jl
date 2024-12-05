"""
    SyntheticFields

Provides implementations of synthetic vector fields.

These can be used to represent a "normal" fluid velocity field which influences the motion
of vortex lines.
"""
module SyntheticFields

export
    SyntheticVectorField,
    UniformVectorField,
    FourierSyntheticVectorField,
    FourierBandVectorField

using StaticArrays: SVector
using Random: Random, AbstractRNG
using LinearAlgebra: ⋅, ×

"""
    SyntheticVectorField{T, N}

Abstract type representing a synthetic vector field in ``N`` dimensions.

Here `T <: AbstractFloat` is the type of the returned values when evaluating the field at a
position.

A field can be evaluated using the `f(x⃗)` syntax, where `f` is a `SyntheticVectorField` and
`x⃗` is a physical location, returning an `SVector{N, T}`. Here `x⃗` can be an `N`-element
tuple or `SVector`.
"""
abstract type SyntheticVectorField{T <: AbstractFloat, N} end

include("uniform.jl")
include("fourier.jl")

end
