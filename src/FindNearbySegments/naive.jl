export NaiveSegmentFinder

"""
    NaiveSegmentFinder <: NearbySegmentFinder
    NaiveSegmentFinder(fs::AbstractVector{<:AbstractFilament})

Initialise nearby segment finder based on a naive searching method.

When using this segment finder, [`nearby_segments`](@ref) iterates through *all* segments in
all filaments `fs`, performing no filtering at all.
One should manually use [`segment_is_close`](@ref) to filter the returned segments while
iterating.
"""
struct NaiveSegmentFinder{
        Filaments <: AbstractVector{<:AbstractFilament},
    } <: NearbySegmentFinder
    fs    :: Filaments
    function NaiveSegmentFinder(fs::AbstractVector{<:AbstractFilament})
        new{typeof(fs)}(copy(fs))
    end
end

function set_filaments!(c::NaiveSegmentFinder, fs)
    resize!(c.fs, length(fs))
    for i ∈ eachindex(c.fs)
        c.fs[i] = fs[i]  # copies filament references ("pointers")
    end
    c
end

struct NaiveSegmentIterator{
        Filament <: AbstractFilament,
        FilamentVector <: AbstractVector{Filament},
    }
    fs :: FilamentVector
end

Base.IteratorSize(::Type{<:NaiveSegmentIterator}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:NaiveSegmentIterator}) = Base.HasEltype()
function Base.eltype(::Type{<:NaiveSegmentIterator{Filament}}) where {Filament}
    S = Segment{Filament}
    @assert isconcretetype(S)
    Tuple{Int, S}  # (filament index, segment)
end

function Base.iterate(it::NaiveSegmentIterator)
    (; fs,) = it
    isempty(fs) && return nothing
    i = firstindex(fs)
    @inbounds f = fs[i]
    j = firstindex(segments(f))
    iterate(it, (i, j))
end

function Base.iterate(it::NaiveSegmentIterator, state::Tuple)
    (; fs,) = it
    i, j = state  # current filament index, current segment index
    @inbounds segs = segments(fs[i])
    while j > lastindex(segs)  # we're done iterating over segments of filament i
        i += 1  # jump to next filament
        if i > lastindex(fs)  # iterated over all filaments
            return nothing
        end
        @inbounds segs = segments(fs[i])
        j = firstindex(segs)
    end
    @inbounds element = (i, segs[j])
    state = (i, j + 1)
    element, state
end

# The iterator should return a tuple (i, segment) where `i` is the index of the filament
# containing the segment (that is, the segment belongs to the filament fs[i]).
nearby_segments(c::NaiveSegmentFinder, x⃗::Vec3) = NaiveSegmentIterator(c.fs)
