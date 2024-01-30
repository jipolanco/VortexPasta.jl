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

function nearby_segments(c::NaiveSegmentFinder, x⃗::Vec3)
    (; fs,) = c
    Filament = eltype(fs)
    @assert Filament <: AbstractFilament
    S = Segment{Filament}
    @assert isconcretetype(S)
    iters = Iterators.map(segments, fs)
    it = Iterators.flatten(iters)
    # @show eltype(it)  # this gives Any, but luckily that doesn't seem to affect type stability
    it
end
