"""
    Diagnostics

Contains tools for computing different diagnostics (total energy, energy spectra, ...) from
simulation data.
"""
module Diagnostics

using ..PaddedArrays: PaddedVector
using ..Filaments:
    Filaments,
    AbstractFilament, ClosedFilament,
    Derivative, UnitTangent, Vec3,
    knots, segments, integrate,
    number_type

using ..BiotSavart: BiotSavart, BiotSavartCache, LongRangeCache, Infinity, ∞

using Bumper: Bumper, @no_escape, @alloc
using LinearAlgebra: ⋅, ×

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const SingleFilamentData = AbstractVector{<:Vec3}
const SetOfFilamentsData = AbstractVector{<:SingleFilamentData}

# Trait indicating whether we can interpolate a list of values, e.g. a vector of velocities
# on filament discretisation points. A filament (or a filament-like object) can be
# interpolated, while a regular vector cannot.
struct IsInterpolable{B}
    IsInterpolable(b::Bool) = new{b}()
end

isinterpolable(::Type{<:AbstractVector}) = IsInterpolable(false)
isinterpolable(::Type{<:AbstractFilament}) = IsInterpolable(true)
isinterpolable(u::AbstractVector) = isinterpolable(typeof(u))  # note: this also applies to filaments (since AbstractFilament <: AbstractVector)

include("energy.jl")
include("helicity.jl")
include("filament_length.jl")
include("vortex_impulse.jl")
include("stretching.jl")
include("spectra.jl")

end
