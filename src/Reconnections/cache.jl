using ..FindNearbySegments:
    FindNearbySegments,
    NearbySegmentFinder,
    NaiveSegmentFinder,
    CellListSegmentFinder,
    set_filaments!,
    segment_is_close,
    nearby_segments
using ..CellLists: CellLists
using OhMyThreads: OhMyThreads

abstract type AbstractReconnectionCache end

distance(c::AbstractReconnectionCache) = distance(criterion(c))

struct ReconnectionCache{
        Criterion <: ReconnectionCriterion,
        Finder <: NearbySegmentFinder,
        ReconnectionInfo,
        Periods <: Tuple{Vararg{Real}},
    } <: AbstractReconnectionCache
    crit     :: Criterion
    finder   :: Finder
    to_reconnect :: Vector{ReconnectionInfo}  # list of segment pairs to be reconnected
    Ls       :: Periods
end

criterion(c::ReconnectionCache) = c.crit
periods(c::ReconnectionCache) = c.Ls
max_passes(c::ReconnectionCache) = max_passes(c.crit)

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

In the second case, the detection of reconnection pairs can follow two different
strategies:

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
    d_reconnect = distance(crit)
    r_cut = 2 * d_reconnect
    finder = if has_nonperiodic_directions
        NaiveSegmentFinder(fs)
    else
        Lmin = min(Ls...)
        M = 2  # number of subdivisions
        r_cut_max = CellLists.max_cutoff_distance(M, Lmin)
        r_cut > r_cut_max && error(
            lazy"""reconnection distance is too large compared to the domain size: d_reconnect / L_min = $(d_reconnect / Lmin).
            Try a smaller reconnection distance or a larger domain."""
        )
        CellListSegmentFinder(fs, r_cut, Ls; nsubdiv = Val(M))
    end
    to_reconnect = let segs = segments(first(fs))
        # Determine ReconnectionInfo type by creating a single candidate
        candidate = ReconnectionCandidate(segs[1], segs[2], firstindex(fs), firstindex(fs))
        info = should_reconnect(crit, candidate; periods = Ls, can_return_nothing = Val(false))::NamedTuple
        reconnect_info = (; candidate, info,)
        ReconnectionInfo = typeof(reconnect_info)
        @assert isconcretetype(ReconnectionInfo)
        Vector{ReconnectionInfo}()
    end
    ReconnectionCache(crit, finder, to_reconnect, Ls)
end

# Used when reconnections are disabled.
struct NullReconnectionCache <: AbstractReconnectionCache end
_init_cache(::NoReconnections, args...) = NullReconnectionCache()
criterion(::NullReconnectionCache) = NoReconnections()

function find_reconnection_pairs!(
        cache::ReconnectionCache,
        fs::AbstractVector{<:AbstractFilament},
        vs::Union{Nothing, AbstractVector{<:AbstractFilament}} = nothing;
        to = TimerOutput(),
    )
    (; finder, to_reconnect,) = cache
    crit = criterion(cache)
    Ls = periods(cache)
    r_cut = distance(cache)
    if crit.use_velocity && vs === nothing
        error("`use_velocity` was set to `true` in the reconnection criterion, but velocity information was not passed to `reconnect!`")
    end
    empty!(to_reconnect)
    T = typeof(r_cut)
    r_crit = T(1.5) * r_cut
    r²_crit = r_crit * r_crit
    Lhs = map(L -> L / 2, Ls)
    @timeit to "set_filaments!" set_filaments!(finder, fs)  # this is needed in particular to initialise cell lists
    lck = ReentrantLock()
    Nf = length(fs)
    # TODO: preallocate arrays?
    ntasks = Threads.nthreads()
    invalidated_filaments = falses(Nf)  # this allows to add filaments at most once (ideally; in parallel that's not always the case, but that's ok for now)
    invalidated_per_thread = [falses(Nf) for _ in 1:ntasks]
    @timeit to "find segment pairs" for (i, f) ∈ pairs(fs)
        # Synchronise per-thread values onto invalidated_filaments.
        for invalidated_local in invalidated_per_thread
            @. invalidated_filaments = invalidated_filaments | invalidated_local
        end
        invalidated_filaments[i] && continue  # skip if this filament has already been invalidated
        @sync for (itask, segments_task) in enumerate(OhMyThreads.index_chunks(eachindex(segments(f)); n = ntasks))
            OhMyThreads.@spawn for i_seg in segments_task
                invalidated_local = invalidated_per_thread[itask]
                copyto!(invalidated_local, invalidated_filaments)  # copy most recent invalidated filaments
                seg_a = Segment(f, i_seg)
                # Since we only compare the *midpoint* of this segment to the extrema of other
                # segments, we add δ/2 (half the segment length) to the critical distance to take
                # into account the case where the point of minimum distance is at the extrema of
                # this segment (and not near the midpoint).
                x⃗ = Filaments.midpoint(seg_a)
                δ² = sum(abs2, f[seg_a.i + 1] - f[seg_a.i])  # ≈ squared segment length
                d²_crit = r²_crit + δ² / 4
                d_crit = @fastmath sqrt(d²_crit)
                for (j, seg_b) ∈ nearby_segments(finder, x⃗)
                    invalidated_local[j] && continue  # skip if filament j has been invalidated
                    # 1. Apply slightly finer filters to determine whether we keep this candidate.
                    # TODO combine these two criteria?
                    @inline segment_is_close(seg_b, x⃗, d_crit, d²_crit, Ls, Lhs) || continue
                    keep_segment_pair(seg_a, seg_b) || continue
                    candidate = ReconnectionCandidate(seg_a, seg_b, i, j)
                    # 2. Now check if the chosen criterion is verified
                    info = should_reconnect(crit, candidate; periods = Ls)
                    info === nothing && continue
                    # 3. Find possibly better candidates
                    # TODO: make sure there are no repetitions...
                    info, candidate = find_better_candidates(info, candidate) do other_candidate
                        should_reconnect(crit, other_candidate; periods = Ls)
                    end
                    (; a, b, filament_idx_a, filament_idx_b,) = candidate
                    @assert info !== nothing
                    reconnect_info = (; candidate, info,)
                    # 4. Check if the candidate satisfies the velocity criterion
                    (; d⃗,) = info  # d⃗ = x⃗ - (y⃗ - p⃗)
                    if crit.use_velocity
                        # Use velocity information to discard some reconnections
                        @assert vs !== nothing   # checked earlier
                        @assert eltype(vs) <: AbstractFilament  # this means it's interpolable
                        v⃗_a = vs[filament_idx_a](a.i, info.ζx)  # assume velocity is interpolable
                        v⃗_b = vs[filament_idx_b](b.i, info.ζy)  # assume velocity is interpolable
                        v_d = d⃗ ⋅ (v⃗_a - v⃗_b)  # separation velocity (should be divided by |d⃗| = sqrt(d²), but we only care about the sign)
                        if v_d > 0  # they're getting away from each other
                            continue  # don't reconnect them
                        end
                    end
                    # 5. If the candidate passed all the tests, we add it to the reconnection list (unless it's already there)
                    invalidated_local[filament_idx_a] = true
                    invalidated_local[filament_idx_b] = true
                    lock(lck) do
                        @inline
                        push!(to_reconnect, reconnect_info)
                    end
                end
            end
        end
    end
    # Sort reconnection candidates according to reconnection distance.
    # This gives higher priority to closest segment pairs (which kinda makes sense).
    # Moreover, it ensures that results are independent of iteration order, which
    # (hopefully) allows to have stable results (which don't vary from one run to another)
    # when running with multiple threads.
    if !isempty(to_reconnect)
        # Note: the x⃗ position (one of the reconnection points) is only used for sorting in
        # case multiple candidates have the same distance d². This is highly improbable in
        # practice, but can happen in highly symmetric cases (in tests in particular!).
        sort!(to_reconnect; by = r -> (r.info.d², r.info.x⃗))
    end
    to_reconnect
end

# If we already found a relatively good reconnection candidate, look at the neighbouring
# segments to see if we can further reduce the distance between the filaments.
@inline function find_better_candidates(f::F, info, c::ReconnectionCandidate) where {F <: Function}
    d²_min = info.d²  # this is the minimum distance we've found until now
    while true
        x, y = c.a, c.b
        x′ = choose_neighbouring_segment(x, info.ζx)
        y′ = choose_neighbouring_segment(y, info.ζy)
        if x′ === x && y′ === y
            break  # the proposed segments are the original ones, so we stop here
        end
        # Note: filaments stay the same, so fields filament_idx_* don't change.
        c′ = ReconnectionCandidate(x′, y′, c.filament_idx_a, c.filament_idx_b)
        info′ = f(c′)
        if info′ === nothing || info′.d² ≥ d²_min
            break  # the previous candidate was better, so we stop here
        end
        # We found a better candidate! But we keep looking just in case.
        info = info′
        c = c′
        d²_min = info.d²
    end
    info, c
end

# This is to make sure we don't reconnect nearly neighbouring segments belonging to the same
# filament.
function keep_segment_pair(a::Segment, b::Segment)
    f, i = a.f, a.i
    g, j = b.f, b.i
    if f === g
        dist_max = 2  # disallow reconnection between segment i and segment j = i ± {0, 1, 2}
        N = length(f)
        dist = abs(i - j)
        if dist ≤ dist_max || (N - dist) ≤ dist_max  # note: N - dist is always ≥ 0
            return false
        end
    end
    true
end
