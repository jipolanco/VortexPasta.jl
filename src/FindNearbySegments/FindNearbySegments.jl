"""
    FindNearbySegments

A module for finding pairs of nearby filament segments.

This is a common operation, used for instance for:

- computing short-range Biot-Savart interactions between filament segments;
- detecting vortex reconnection candidates.
"""
module FindNearbySegments

export set_filaments!, nearby_segments, segment_is_close

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

"""
    segment_is_close(s::Segment, x⃗, r_cut::Real, r²_cut::Real, Ls, Lhs) -> Bool

Determine whether segment `s` is close to the point `x⃗`.

Here `r_cut` is the critical distance below which `s` and `x⃗` will be considered to be
"close", and `r²_cut = r_cut^2` is its squared value.

Moreover, `Ls = (Lx, Ly, Lz)` contains the domain period along each dimension (components
can be `Infinity()` for infinite non-periodic domains), and `Lhs = (Lx/2, Ly/2, Lz/2)`
contains half the domain periods.
"""
function segment_is_close(s::Segment, x⃗, r_cut::Real, r²_cut::Real, Ls, Lhs)
    r⃗_a = let  # separation with first node of the segment
        y⃗ = @inbounds s.f[s.i]
        deperiodise_separation(y⃗ - x⃗, Ls, Lhs)
    end
    r⃗_b = let  # separation with second node of the segment
        y⃗ = @inbounds s.f[s.i + 1]
        deperiodise_separation(y⃗ - x⃗, Ls, Lhs)
    end

    # Skip segment if both extrema are too far from the point of interest.
    if any(>(r_cut), r⃗_a) && any(>(r_cut), r⃗_b)
        return false
    end

    # Second (more precise) filter: look at the actual distances.
    if sum(abs2, r⃗_a) > r²_cut && sum(abs2, r⃗_b) > r²_cut
        return false
    end

    # The current segment is sufficiently close to x⃗, so we return it.
    true
end

include("naive.jl")
include("cell_lists.jl")

end
