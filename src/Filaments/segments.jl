"""
    Segments{Filament <: AbstractFilament}

Convenience type allowing to iterate over the segments of a filament.
"""
struct Segments{Filament <: AbstractFilament}
    f :: Filament
end

"""
    segments(f::AbstractFilament) -> Segments(f)

Create a [`Segments`](@ref) object for iterating over the segments of a filament.
"""
@inline segments(f::AbstractFilament) = Segments(f)

"""
    Base.length(s::Segments{<:AbstractFilament}) -> Int

Return the number of segments in a filament.
"""
Base.length(s::Segments) = _length(s.f, s)
_length(f::ClosedFilament, ::Segments) = length(f)

"""
    Base.eachindex(s::Segments{<:AbstractFilament}) -> AbstractRange

Return the indices associated to the segments of a filament.
"""
@inline Base.eachindex(s::Segments) = _eachindex(s.f, s)
_eachindex(f::ClosedFilament, ::Segments) = eachindex(f)

@inline Base.firstindex(s::Segments) = first(eachindex(s))
@inline Base.lastindex(s::Segments) = last(eachindex(s))
