"""
    FilamentIO

Module defining tools for reading and writing filament data.
"""
module FilamentIO

using ..PaddedArrays: PaddedVector, pad_periodic!
using ..Filaments
using ..Filaments: DiscretisationMethod, GeometricQuantity, Vec3

# Return i-th on-segment location ζ ∈ [0, 1] for a given level of refinement.
# Here i ∈ 1:refinement.
on_segment_location(i, refinement) = (i - 1) / refinement

include("json_vtk_series.jl")
include("hdf5.jl")
    
end
