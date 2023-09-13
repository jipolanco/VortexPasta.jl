"""
    Diagnostics

Contains tools for computing different diagnostics (total energy, energy spectra, ...) from
simulation data.
"""
module Diagnostics

using ..Filaments: Filaments, Derivative, Vec3, knots, segments, integrate
using LinearAlgebra: ⋅

const SingleFilamentData = AbstractVector{<:Vec3}
const SetOfFilamentsData = AbstractVector{<:SingleFilamentData}

include("energy.jl")

end
