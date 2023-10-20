export NaiveSegmentFinder

struct NaiveSegmentIterator{
        S <: Segment,
        Filaments <: AbstractVector{<:AbstractFilament},
        CutoffDistance <: Real,
        Periods <: NTuple{3, Real},
        HalfPeriods <: NTuple{3, Real},
        Position <: Vec3{<:Real},
    } <: NearbySegmentIterator{S}
    fs     :: Filaments
    r_cut  :: CutoffDistance
    r²_cut :: CutoffDistance
    Ls     :: Periods
    Lhs    :: HalfPeriods
    x⃗      :: Position

    function NaiveSegmentIterator{S}(fs, r_cut, r²_cut, Ls, Lhs, x⃗) where {S}
        new{S, typeof(fs), typeof(r_cut), typeof(Ls), typeof(Lhs), typeof(x⃗)}(
            fs, r_cut, r²_cut, Ls, Lhs, x⃗,
        )
    end
end

@inline function _initial_state(it::NaiveSegmentIterator)
    (; fs, x⃗, Ls, Lhs, r_cut,) = it
    n = firstindex(fs)
    f = first(fs)
    i = firstindex(segments(f)) - 1
    r⃗_b = let  # separation vector to first node of first filament
        y⃗ = f[i + 1]
        deperiodise_separation(y⃗ - x⃗, Ls, Lhs)
    end
    # We discard a pair interaction if the distance *along a Cartesian direction* is larger
    # than r_cut. This is less precise (but cheaper than) using the actual squared distance.
    is_outside_range_b = any(>(r_cut), r⃗_b)
    (; n, i, r⃗_b, is_outside_range_b,)
end

# We mark this function with :terminates_locally just to ensure the compiler that the
# function always terminates (despite the `while true`).
@inline Base.@assume_effects :terminates_locally function Base.iterate(
        it::NaiveSegmentIterator,
        state = _initial_state(it),
    )
    (; fs, r_cut, r²_cut, Ls, Lhs, x⃗,) = it
    (; n, i, r⃗_b, is_outside_range_b,) = state
    State = typeof(state)

    f = @inbounds fs[n]

    # Advance until we find a filament segment which is within the cutoff radius (or until
    # we have explored all filament segments).
    while true
        if i == lastindex(segments(f))
            n += 1
            if n == lastindex(fs) + 1
                return nothing  # we iterated over all filaments, so we're done
            end
            f = @inbounds fs[n]  # next filament
            i = firstindex(segments(f))
        else
            i += 1
        end

        r⃗_a = r⃗_b
        is_outside_range_a = is_outside_range_b

        r⃗_b = let  # separation vector to first node of first filament
            y⃗ = @inbounds f[i + 1]  # the current segment is f[i]..f[i + 1]
            deperiodise_separation(y⃗ - x⃗, Ls, Lhs)
        end
        is_outside_range_b = any(>(r_cut), r⃗_b)

        # Skip this segment if its two limits are too far from x⃗.
        if is_outside_range_a && is_outside_range_b
            continue
        end

        # Second (more precise) filter: look at the actual distances.
        if sum(abs2, r⃗_a) > r²_cut && sum(abs2, r⃗_b) > r²_cut
            continue
        end

        # The current segment is sufficiently close to x⃗, so we return it.
        ret = Segment(f, i) :: eltype(it)  # check that we return the type promised by eltype
        state = (; n, i, r⃗_b, is_outside_range_b,) :: State
        return ret, state
    end
end

## ================================================================================ ##

"""
    NaiveSegmentFinder <: NearbySegmentFinder
    NaiveSegmentFinder(
        fs::AbstractVector{<:AbstractFilament},
        r_cut::Real,
        Ls::NTuple{3, Real},
    )

Initialise nearby segment finder based on a naive searching method.

For a given position of interest ``\\bm{x}``, it computes the distance between this position
and *all* filament segments, returning those segments which are closer than the cut-off
distance `r_cut`.
This can become inefficient when there are many filaments.

# Mandatory arguments

- `fs`: vector of filaments;
- `r_cut::Real`: cut-off distance;
- `Ls::NTuple{3, Real}`: domain period in each direction.
  Can be `Ls = (Infinity(), Infinity(), Infinity())` for an infinite non-periodic domain.
"""
struct NaiveSegmentFinder{
        Filaments <: AbstractVector{<:AbstractFilament},
        CutoffDistance <: Real,
        Periods <: Tuple{Vararg{Real}},
    } <: NearbySegmentFinder
    fs    :: Filaments
    r_cut :: CutoffDistance
    Ls    :: Periods

    function NaiveSegmentFinder(fs::AbstractVector{<:AbstractFilament}, r_cut, Ls)
        new{typeof(fs), typeof(r_cut), typeof(Ls)}(copy(fs), r_cut, Ls)
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
    (; fs, r_cut, Ls,) = c
    Lhs = map(L -> L / 2, Ls)  # half periods
    Filament = eltype(fs)
    @assert Filament <: AbstractFilament
    S = Segment{Filament}
    @assert isconcretetype(S)
    NaiveSegmentIterator{S}(fs, r_cut, r_cut^2, Ls, Lhs, x⃗)
end
