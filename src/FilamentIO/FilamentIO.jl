"""
    FilamentIO

Module defining tools for reading and writing filament data.
"""
module FilamentIO

using ..PaddedArrays: PaddedVector, pad_periodic!
using ..Filaments
using ..Filaments: DiscretisationMethod, Vec3

include("hdf5.jl")
    
end
