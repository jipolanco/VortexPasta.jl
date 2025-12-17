using ..FindNearbySegments:
    FindNearbySegments,
    NearbySegmentFinder,
    NaiveSegmentFinder,
    CellListSegmentFinder,
    set_filaments!,
    segment_is_close,
    nearby_segments

using Accessors: @set

"""
    ReconnectBasedOnDistance <: ReconnectionCriterion
    ReconnectBasedOnDistance(d_crit; max_passes = 1, use_velocity = (max_passes == 1), decrease_length = true, cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

# Optional keyword arguments

- `use_velocity`: if `true` (default), use filament velocity information in reconnections.
  For now, this is used to discard a reconnection between two points if they are
  instantaneously getting away from each other. If `max_passes > 1`, then `use_velocity` **must**
  be set to `false`.

- `max_passes`: maximum number of scans through all segment pairs to detect and perform
  reconnections. In a single pass, a filament can only be reconnected at most once.
  Therefore, this parameter can be useful to make sure that all required reconnections have
  been performed. **Setting `max_passes > 1` also requires `use_velocity = false`.**

- `decrease_length`: if `true` (default), a reconnection will only be performed
  if it will decrease the total filament length. Since, for vortices, the total
  energy is roughly related to the vortex length, this means that reconnections
  should always tend to dissipate energy.

- `cos_max`: allows to disable reconnections of nearly parallel segments. Two segments
  are considered to be "nearly parallel" if `cos(θ) > cos_max`.
  The default value `cos_max = 0.97` disables reconnections when the angle between lines
  is ``θ < \\arccos(0.97) ≈ 14°``.
  Note that the angle ``θ`` is signed (it takes values in ``[-1, 1]``).
  Negative angles mean that the segments are antiparallel, and in this case reconnections are always performed.
"""
struct ReconnectBasedOnDistance <: ReconnectionCriterion
    dist       :: Float64
    dist_sq    :: Float64
    cos_max    :: Float64
    cos_max_sq :: Float64
    max_passes :: Int
    use_velocity    :: Bool
    decrease_length :: Bool
    function ReconnectBasedOnDistance(dist; cos_max = 0.97, max_passes = 1, use_velocity = (max_passes == 1), decrease_length = true)
        max_passes < 1 && throw(ArgumentError("ReconnectBasedOnDistance: max_passes should be >= 1"))
        max_passes > 1 && use_velocity && throw(ArgumentError("ReconnectBasedOnDistance: use_velocity can only be activated if max_passes = 1"))
        new(dist, dist^2, cos_max, cos_max^2, max_passes, use_velocity, decrease_length)
    end
end

distance(c::ReconnectBasedOnDistance) = c.dist
max_passes(c::ReconnectBasedOnDistance) = c.max_passes
require_interpolated_velocity(c::ReconnectBasedOnDistance) = c.use_velocity

function Base.show(io::IO, c::ReconnectBasedOnDistance)
    print(io, "ReconnectBasedOnDistance($(c.dist); cos_max = $(c.cos_max), max_passes = $(c.max_passes), use_velocity = $(c.use_velocity), decrease_length = $(c.decrease_length))")
end

"""
    should_reconnect(
        c::ReconnectBasedOnDistance,
        fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
    ) -> Union{Nothing, NamedTuple}

Check whether two filaments should reconnect according to the chosen criterion.

Checks for a possible reconnection between filament segments `fx[i:i+1]` and `fy[j:j+1]`.

If the filament segments should reconnect, this function returns a `NamedTuple`
with reconnection information, which includes in particular all the fields
returned by [`find_min_distance`](@ref).

Otherwise, returns `nothing` if the filament segments should not reconnect.
"""
function should_reconnect end

should_reconnect(c::ReconnectBasedOnDistance, candidate::ReconnectionCandidate; kws...) =
    should_reconnect(c, candidate.a, candidate.b; kws...)
should_reconnect(c::ReconnectBasedOnDistance, a::Segment, b::Segment; kws...) =
    should_reconnect(c, a.f, b.f, a.i, b.i; kws...)

function should_reconnect(
        c::ReconnectBasedOnDistance, fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
        can_return_nothing = Val(true),  # this is set to false to obtain the return type when constructing a ReconnectBasedOnDistanceCache
    )
    (; dist_sq, cos_max_sq, decrease_length,) = c

    min_dist = find_min_distance(fx, fy, i, j; periods)
    (; d⃗, p⃗, ζx, ζy,) = min_dist
    d² = sum(abs2, d⃗)
    if can_return_nothing === Val(true)
        d² > dist_sq && return nothing  # don't reconnect
    end

    # Make sure that reconnections reduce the total length (makes sense energetically for vortices).
    length_before = norm(fx[i + 1] - fx[i]) + norm(fy[j + 1] - fy[j])
    length_after = norm(fy[j + 1] - fx[i] - p⃗) + norm(fx[i + 1] - fy[j] + p⃗)
    if decrease_length && can_return_nothing === Val(true)
        length_after > length_before && return nothing
    end

    X′ = fx(i, ζx, Derivative(1))
    Y′ = fy(j, ζy, Derivative(1))

    # Return the output of find_min_distance + other stuff if segments should reconnect.
    info = (;
        min_dist...,
        d², length_before, length_after,
    )

    xy = X′ ⋅ Y′
    xy < 0 && return info  # always reconnect antiparallel vortices

    cos² = (xy * xy) / (sum(abs2, X′) * sum(abs2, Y′))

    if cos² < cos_max_sq || can_return_nothing !== Val(true)
        info
    else
        nothing
    end
end

struct ReconnectBasedOnDistanceCache{
        Criterion <: ReconnectBasedOnDistance,
        Finder <: NearbySegmentFinder,
        ReconnectionInfo,
        Periods <: Tuple{Vararg{Real}},
    } <: AbstractReconnectionCache
    crit     :: Criterion
    finder   :: Finder
    to_reconnect :: Vector{ReconnectionInfo}  # list of segment pairs to be reconnected
    Ls       :: Periods
end

criterion(c::ReconnectBasedOnDistanceCache) = c.crit
periods(c::ReconnectBasedOnDistanceCache) = c.Ls

function _init_cache(crit::ReconnectBasedOnDistance, fs, Ls)
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
    ReconnectBasedOnDistanceCache(crit, finder, to_reconnect, Ls)
end

function find_reconnection_pairs!(
        cache::ReconnectBasedOnDistanceCache,
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

function _reconnect_pass!(callback::F, ret_base, cache::ReconnectBasedOnDistanceCache, fs, vs; to) where {F <: Function}
    @timeit to "find reconnection pairs" begin
        to_reconnect = find_reconnection_pairs!(cache, fs, vs; to)
    end
    if isempty(to_reconnect)
        ret = @set ret_base.npasses += 1
        return ret
    end
    (; reconnection_count, reconnection_length_loss, filaments_removed_count, filaments_removed_length,) = ret_base
    Nf = length(fs)
    invalidated_filaments = falses(Nf)  # if invalidated_filaments[i] is true, it means fs[i] can no longer be reconnected
    @timeit to "reconnect pairs" for reconnect_info ∈ to_reconnect
        (; candidate, info,) = reconnect_info
        (; a, b, filament_idx_a, filament_idx_b,) = candidate
        @timeit to "reconnect" if filament_idx_a === filament_idx_b
            if invalidated_filaments[filament_idx_a]
                continue  # don't reconnect this filament if it was already reconnected earlier
            end
            # Reconnect filament with itself => split filament into two
            @assert a.i ≠ b.i
            nremoved, Lremoved = reconnect_with_itself!(callback, fs, a.f, a.i, b.i, info)
            filaments_removed_count += nremoved
            filaments_removed_length += Lremoved
            invalidated_filaments[filament_idx_a] = true
        else
            # Reconnect two different filaments => merge them into one
            if invalidated_filaments[filament_idx_a] || invalidated_filaments[filament_idx_b]
                continue  # don't reconnect if one of these filaments was already reconnected earlier
            end
            reconnect_with_other!(callback, fs, a.f, b.f, a.i, b.i, info)
            invalidated_filaments[filament_idx_a] = true
            invalidated_filaments[filament_idx_b] = true
        end
        reconnection_length_loss += info.length_before - info.length_after
        reconnection_count += 1
    end
    (;
        reconnection_count,
        reconnection_length_loss,
        filaments_removed_count,
        filaments_removed_length,
        npasses = ret_base.npasses + 1,
    ) :: typeof(ret_base)
end

# Here ζ ∈ [0, 1] is the location within the segment where the minimum distance was found by
# `should_reconnect`. If ζ = 0, the minimum location was found at the segment starting
# point, which means that we might find a smaller distance if we look at the previous
# segment. Similar thing if ζ = 1.
@inline function choose_neighbouring_segment(s::Segment, ζ)
    (; f, i,) = s
    if ζ == 0
        j = ifelse(i == firstindex(f), lastindex(f), i - 1)  # basically j = i - 1 (except when i = 1)
        Segment(f, j)
    elseif ζ == 1
        j = ifelse(i == lastindex(f), firstindex(f), i + 1)
        Segment(f, j)
    else
        s
    end
end
