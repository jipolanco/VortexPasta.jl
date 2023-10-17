"""
    FindNearbySegments

A module for finding pairs of nearby filament segments.

This is a common operation, used for instance for:

- computing short-range Biot-Savart interactions between filament segments;
- detecting vortex reconnection candidates.
"""
module FindNearbySegments

using ..Filaments:
    Filaments,
    AbstractFilament,
    Segment,
    Vec3,
    segments,
    deperiodise_separation

# This may be overloaded by different backends (but it's not necessary).
abstract type NearbySegmentIterator{S <: Segment} end

# These are needed e.g. by collect(it::NearbySegmentIterator)
Base.IteratorSize(::Type{<:NearbySegmentIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:NearbySegmentIterator{S}}) where {S <: Segment} = S

abstract type NearbySegmentFinder end

"""
    set_filaments!(c::NearbySegmentFinder, fs::AbstractVector{<:AbstractFilament})

Store (and optionally process) the list of filaments.

This must be called to set the filament list before using [`nearby_segments`](@ref).
"""
function set_filaments! end

"""
    nearby_segments(c::NearbySegmentFinder, x⃗::Vec3)

Return an iterator over the segments that are "close" to the location `x⃗`.

A segment is considered to be close to `x⃗` if the minimum[^mindist] distance between `x⃗`
and the segment midpoint is smaller than the cutoff distance ``r_{\\text{cut}}``.

Typical usage:

```julia
x⃗ = Vec3(0.1, 0.3, 0.2)
for segment ∈ nearby_segments(c, x⃗)
    # Get the filament `f` and the index `i` of the segment within the filament.
    # The segment is between the filament nodes `f[i]` and `f[i + 1]`.
    (; f, i,) = segment
    # Do something with the segment...
end
```

[^mindist]: When periodicity is enabled, the relevant distance is the *minimum* distance
            between the two points, after considering all their periodic images.
"""
function nearby_segments end

include("naive.jl")
include("cell_lists.jl")

end
