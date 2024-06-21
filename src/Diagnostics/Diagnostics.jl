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

include("energy.jl")
include("helicity.jl")
include("filament_length.jl")
include("vortex_impulse.jl")
include("stretching.jl")
include("spectra.jl")

end
