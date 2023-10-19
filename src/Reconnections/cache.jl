using ..FindNearbySegments:
    FindNearbySegments,
    NearbySegmentFinder,
    NaiveSegmentFinder,
    CellListSegmentFinder,
    set_filaments!,
    nearby_segments
using ..CellLists: CellLists  # for docs only

abstract type AbstractReconnectionCache end

distance(c::AbstractReconnectionCache) = distance(criterion(c))

struct ReconnectionCache{
        Criterion <: ReconnectionCriterion,
        Finder <: NearbySegmentFinder,
        Candidate <: ReconnectionCandidate,
        Periods <: Tuple{Vararg{Real}},
    } <: AbstractReconnectionCache
    crit       :: Criterion
    finder     :: Finder
    candidates :: Vector{Union{Nothing, Candidate}}
    Ls         :: Periods
end

criterion(c::ReconnectionCache) = c.crit
periods(c::ReconnectionCache) = c.Ls

function invalidate_candidates!(cache::ReconnectionCache, flist::Vararg{AbstractFilament})
    (; candidates,) = cache
    for (i, candidate) ∈ pairs(candidates)
        candidate === nothing && continue
        (; a, b,) = candidate
        for f ∈ flist
            if f === a.f || f === b.f
                candidates[i] = nothing  # invalidate candidate
            end
        end
    end
    candidates
end

"""
    Reconnections.init_cache(
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament},
        Ls::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) -> AbstractReconnectionCache

Initialise reconnection cache.

Required arguments are a reconnection criterion `crit` and the domain dimensions (or *periods*) `Ls`.

The type of cache will vary depending on the inputs:

- if `crit = NoReconnections()`, then reconnections are disabled and this returns a
  `NullReconnectionCache`;
- otherwise, it returns a `ReconnectionCache`.

In the second case, the detection of reconnection candidates can follow two different strategies:

- if the domain is infinite (default), a naive iteration across all filament segment pairs
is performed;
- if the domain is periodic, an efficient cell lists algorithm is used (see the [`CellLists`](@ref) module).
"""
function init_cache(
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament},
        Ls::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )
    _init_cache(crit, fs, Ls)
end

function _init_cache(crit::ReconnectionCriterion, fs, Ls)
    has_nonperiodic_directions = any(L -> L === Infinity(), Ls)
    # Note: we make the cutoff distance larger than the actual critical distance, since this
    # distance is only used to compare the segment *midpoints*.
    r_cut = 2 * distance(crit)
    finder = if has_nonperiodic_directions
        NaiveSegmentFinder(fs, r_cut, Ls)
    else
        CellListSegmentFinder(fs, r_cut, Ls; nsubdiv = Val(2))
    end
    candidates = let
        FilamentType = eltype(fs)
        S = Segment{FilamentType}
        T = ReconnectionCandidate{S}
        @assert isconcretetype(T)
        Union{T, Nothing}[]
    end
    ReconnectionCache(crit, finder, candidates, Ls)
end

# Used when reconnections are disabled.
struct NullReconnectionCache <: AbstractReconnectionCache end
_init_cache(::NoReconnections, args...) = NullReconnectionCache()
criterion(::NullReconnectionCache) = NoReconnections()

function find_reconnection_candidates!(
        cache::ReconnectionCache,
        fs::AbstractVector{<:AbstractFilament};
        to = TimerOutput(),
    )
    (; finder, candidates,) = cache
    empty!(candidates)
    r_cut = distance(cache)
    r²_crit = 4 * r_cut * r_cut
    Ls = periods(cache)
    Lhs = map(L -> L / 2, Ls)
    function check_distance(r⃗_in)
        r⃗ = deperiodise_separation(r⃗_in, Ls, Lhs)
        r² = sum(abs2, r⃗)
        r² < r²_crit
    end
    @timeit to "set_filaments!" set_filaments!(finder, fs)  # this is needed in particular to initialise cell lists
    @timeit to "add candidates" for f ∈ fs, seg_a ∈ segments(f)
        x⃗ = Filaments.midpoint(seg_a)
        for seg_b ∈ nearby_segments(finder, x⃗)
            # Slightly finer filter to determine whether we keep this candidate.
            keep_segment_pair(check_distance, seg_a, seg_b) || continue
            push!(candidates, ReconnectionCandidate(seg_a, seg_b))
        end
    end
    candidates
end

# This is to make sure we don't reconnect nearly neighbouring segments belonging to the same
# filament.
function keep_segment_pair(check_distance::F, a::Segment, b::Segment) where {F <: Function}
    f, i = a.f, a.i
    g, j = b.f, b.i
    if f === g
        i === j && return false  # same segment
        dist_max = 2  # disallow reconnection between segment i and segment j = i ± {0, 1, 2}
        N = length(f)
        Nh = N >> 1
        i′, j′ = ifelse(i < j, (i, j), (j, i))
        while j′ - i′ > Nh
            i′, j′ = j′ - N, i′
        end
        if j′ - i′ ≤ dist_max
            return false
        end
    end
    any((
        check_distance(f[i + 0] - g[j + 0]),
        check_distance(f[i + 1] - g[j + 0]),
        check_distance(f[i + 0] - g[j + 1]),
        check_distance(f[i + 1] - g[j + 1]),
    ))
end
