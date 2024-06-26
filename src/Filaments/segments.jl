## SegmentIterator type

"""
    SegmentIterator{Filament <: AbstractFilament}

Convenience type allowing to iterate over the segments of a filament.
"""
struct SegmentIterator{Filament <: AbstractFilament}
    f :: Filament
end

"""
    segments(f::AbstractFilament) -> SegmentIterator(f)

Create a [`SegmentIterator`](@ref) object for iterating over the segments of a filament.
"""
@inline segments(f::AbstractFilament) = SegmentIterator(f)

"""
    Base.length(s::SegmentIterator{<:AbstractFilament}) -> Int

Return the number of segments in a filament.
"""
Base.length(s::SegmentIterator) = _length(s.f, s)
_length(f::ClosedFilament, ::SegmentIterator) = length(f)

"""
    Base.eachindex(s::SegmentIterator{<:AbstractFilament}) -> AbstractRange

Return the indices associated to the segments of a filament.
"""
@inline Base.eachindex(s::SegmentIterator) = _eachindex(s.f, s)
_eachindex(f::ClosedFilament, ::SegmentIterator) = eachindex(f)

@inline Base.firstindex(s::SegmentIterator) = first(eachindex(s))
@inline Base.lastindex(s::SegmentIterator) = last(eachindex(s))

# This is used by pairs(s::SegmentIterator).
@inline Base.keys(s::SegmentIterator) = eachindex(s)

@inline function Base.iterate(it::SegmentIterator, i = firstindex(it))
    if i ≤ lastindex(it)
        @inbounds(it[i]), i + 1
    else
        nothing
    end
end

Base.eltype(::Type{<:SegmentIterator{F}}) where {F} = Segment{F}
Base.eltype(s::SegmentIterator) = eltype(typeof(s))
Base.@propagate_inbounds Base.getindex(s::SegmentIterator, i::Integer) = Segment(s.f, i)


## Segment type

"""
    Segment{<:Filament}
    Segment(f::AbstractFilament, i::Integer)

Represents a single filament segment.

The segment goes from nodes `f[i]` to `f[i + 1]`.
"""
struct Segment{Filament <: AbstractFilament}
    f :: Filament
    i :: Int
end

"""
    (s::Segment)(ζ, [args...])

Evaluate quantity along a segment.

Here ``0 ≤ ζ ≤ 1`` is a location within the segment.

This is the equivalent of `s.f(s.i, ζ, args...)`. As an example, `s(ζ)` simply computes a
coordinate within the segment, while `s(ζ, UnitTangent())` computes the unit tangent at that
coordinate.

See [`AbstractFilament`](@ref) for more details on what can be computed. 
"""
function (s::Segment)(ζ::Number, args...)
    s.f(s.i, ζ, args...)
end

startpoint(s::Segment) = @inbounds s.f[s.i]
endpoint(s::Segment) = @inbounds s.f[s.i + 1]

"""
    midpoint(s::Segment) -> Vec3

Return an estimation of the segment midpoint (prioritising performance over accuracy).
"""
@inline function midpoint(s::Segment)
    # s.f(s.i, 0.5)  # "exact" midpoint
    (startpoint(s) + endpoint(s)) ./ 2  # approximation, usually faster
end

"""
    Filaments.segment_length(s::Segment; quad = nothing)

Estimate length of a filament segment.

One may pass a quadrature rule as `quad` for better accuracy.
Otherwise, if `quad = nothing`, this simply returns the straight distance between the two
segment extremities.
"""
segment_length(s::Segment; quad = nothing) = _segment_length(quad, s)

_segment_length(::Nothing, s::Segment) = norm(startpoint(s) - endpoint(s))

function _segment_length(quad::AbstractQuadrature, s::Segment)
    integrate(s, quad) do s, ζ
        norm(s(ζ, Derivative(1)))
    end
end
