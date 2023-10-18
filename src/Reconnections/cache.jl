using ..FindNearbySegments:
    FindNearbySegments,
    NearbySegmentFinder,
    NaiveSegmentFinder,
    CellListSegmentFinder,
    set_filaments!,
    nearby_segments

abstract type AbstractReconnectionCache end

distance(c::AbstractReconnectionCache) = distance(criterion(c))

struct ReconnectionCandidate{S <: Segment}
    a :: S
    b :: S
end

struct ReconnectionCache{
        Criterion <: ReconnectionCriterion,
        Finder <: NearbySegmentFinder,
        Candidate <: ReconnectionCandidate,
        Periods <: Tuple{Vararg{Real}},
    } <: AbstractReconnectionCache
    crit       :: Criterion
    finder     :: Finder
    candidates :: Vector{Candidate}
    Ls         :: Periods
end

criterion(c::ReconnectionCache) = c.crit
periods(c::ReconnectionCache) = c.Ls

function remove_filaments_from_candidates!(cache::ReconnectionCache, flist::Vararg{AbstractFilament})
    (; candidates,) = cache
    inds = reverse(eachindex(candidates))
    for i ∈ inds
        (; a, b,) = candidates[i]
        if any(g -> g === a.f || g === b.f, flist)
            popat!(candidates, i)
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
    r_cut = 4 * distance(crit)
    finder = if has_nonperiodic_directions
        NaiveSegmentFinder(fs, r_cut, Ls)
    else
        CellListSegmentFinder(fs, r_cut, Ls; nsubdiv = Val(1))
    end
    candidates = let
        FilamentType = eltype(fs)
        S = Segment{FilamentType}
        T = ReconnectionCandidate{S}
        @assert isconcretetype(T)
        T[]
    end
    ReconnectionCache(crit, finder, candidates, Ls)
end

# Used when reconnections are disabled.
struct NullReconnectionCache <: AbstractReconnectionCache end
_init_cache(::NoReconnections, args...) = NullReconnectionCache()
criterion(::NullReconnectionCache) = NoReconnections()

function find_reconnection_candidates!(
        cache::ReconnectionCache,
        fs::AbstractVector{<:AbstractFilament},
    )
    (; finder, candidates,) = cache
    empty!(candidates)
    set_filaments!(finder, fs)  # this is needed in particular to initialise cell lists
    for f ∈ fs, seg_a ∈ segments(f)
        x⃗ = Filaments.midpoint(seg_a)
        for seg_b ∈ nearby_segments(finder, x⃗)
            seg_a === seg_b && continue
            # For now we simply include the pair of segments, but we may want to apply a few
            # more filters to reduce the number of candidates.
            push!(candidates, ReconnectionCandidate(seg_a, seg_b))
        end
    end
    candidates
end
