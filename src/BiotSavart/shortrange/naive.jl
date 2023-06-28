export NaiveShortRangeBackend

"""
    NaiveShortRangeBackend <: ShortRangeBackend

Naive computation of short-range interactions.
"""
struct NaiveShortRangeBackend <: ShortRangeBackend end

struct NaiveShortRangeCache{
        FilamentsRef <: Ref{<:AbstractVector{<:AbstractFilament}},
        Params <: ParamsShortRange,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    fs_ref :: FilamentsRef
    params :: Params
    to     :: Timer
end

function init_cache_short(
        ::ParamsCommon, params::ParamsShortRange{<:NaiveShortRangeBackend},
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    )
    NaiveShortRangeCache(Ref(fs), params, to)
end

function set_filaments!(c::NaiveShortRangeCache, fs)
    c.fs_ref[] = fs
    c
end

struct NaiveSegmentIterator{
        S <: Segment,
        Filaments <: AbstractVector{<:AbstractFilament},
        CutoffRadius <: Real,
        Periods <: NTuple{3, Real},
        HalfPeriods <: NTuple{3, Real},
        Position <: Vec3{<:Real},
    } <: NearbySegmentIterator{S}
    fs     :: Filaments
    r_cut  :: CutoffRadius
    r²_cut :: CutoffRadius
    Ls     :: Periods
    Lhs    :: HalfPeriods
    x⃗      :: Position

    function NaiveSegmentIterator{S}(fs, r_cut, r²_cut, Ls, Lhs, x⃗) where {S}
        new{S, typeof(fs), typeof(r_cut), typeof(Ls), typeof(Lhs), typeof(x⃗)}(
            fs, r_cut, r²_cut, Ls, Lhs, x⃗,
        )
    end
end

function nearby_segments(c::NaiveShortRangeCache, x⃗::Vec3)
    (; params, fs_ref,) = c
    (; common, rcut,) = params
    (; Ls,) = common
    Lhs = map(L -> L / 2, Ls)  # half periods
    Filament = eltype(eltype(fs_ref))
    @assert Filament <: AbstractFilament
    S = Segment{Filament}
    @assert isconcretetype(S)
    NaiveSegmentIterator{S}(fs_ref[], rcut, rcut^2, Ls, Lhs, x⃗)
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
