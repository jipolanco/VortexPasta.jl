"""
    Diagnostics

Contains tools for computing different diagnostics (total energy, energy spectra, ...) from
simulation data.
"""
module Diagnostics

using ..Filaments:
    Filaments,
    AbstractFilament, Derivative, Vec3,
    knots, segments, integrate

using ..BiotSavart: BiotSavartCache, LongRangeCache, Infinity, ∞

using LinearAlgebra: ⋅, ×

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const SingleFilamentData = AbstractVector{<:Vec3}
const SetOfFilamentsData = AbstractVector{<:SingleFilamentData}

include("energy.jl")
include("filament_length.jl")
include("vortex_impulse.jl")
include("spectra.jl")

end
