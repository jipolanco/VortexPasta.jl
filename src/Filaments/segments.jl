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
    midpoint(s::Segment) -> Vec3

Return an estimation of the segment midpoint (prioritising performance over accuracy).
"""
@inline function midpoint(s::Segment)
    (; f, i,) = s
    # f(i, 0.5)  # "exact" midpoint
    @inbounds (f[i] .+ f[i + 1]) ./ 2  # approximation, usually faster
end
