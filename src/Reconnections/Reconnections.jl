"""
    Reconnections

Module for dealing with the reconnection of filaments.
"""
module Reconnections

export ReconnectionCriterion,
       NoReconnections,
       ReconnectBasedOnDistance,
       reconnect!

using ..BasicTypes: Infinity

using ..Filaments: Filaments,
                   AbstractFilament,
                   Segment,
                   Derivative,
                   segments,
                   split!, merge!,
                   update_coefficients!,
                   find_min_distance,
                   filament_length,
                   deperiodise_separation,
                   number_type,
                   check_nodes

using TimerOutputs
using LinearAlgebra: ⋅, norm

struct ReconnectionCandidate{S <: Segment}
    a :: S
    b :: S
end

include("criteria.jl")
include("cache.jl")

"""
    reconnect!(
        [callback::Function],
        cache::AbstractReconnectionCache,
        fs::AbstractVector{<:AbstractFilament},
    ) -> NamedTuple

Perform filament reconnections according to chosen criterion.

Note that, when a filament self-reconnects, this creates new filaments, which
are appended at the end of `fs`.

Moreover, this function will remove reconnected filaments if their number of nodes is too small
(typically ``< 3``, see [`check_nodes`](@ref)).

## Returns

This function returns a `NamedTuple` with fields:

- `reconnection_count`: **number of reconnections**;

- `filament_removed_count`: **number of filaments removed after reconnection**. This happens
  when a filament reconnects with itself, splitting into two new filaments, and one or both
  of these filaments have an insufficient number of discretisation nodes;

- `filament_removed_length`: **length of removed filaments**. We use a straight segment
  approximation (no quadratures) to estimate the filament lengths.

## Callback function

Optionally, one may pass a callback which will be called whenever the vector of
filaments `fs` is modified. Its signature must be the following:

    callback(f::AbstractFilament, i::Int, mode::Symbol)

where `f` is the modified filament, `i` is its index in `fs`, and `mode` is one of:

- `:modified` if the filament `fs[i]` was modified;
- `:appended` if the filament was appended at the end of `fs` (at index `i`);
- `:removed` if the filament previously located at index `i` was removed.

"""
function reconnect!(
        callback::F,
        cache::AbstractReconnectionCache,
        fs::AbstractVector{<:AbstractFilament};
        to::TimerOutput = TimerOutput(),
    ) where {F <: Function}
    T = number_type(fs)
    reconnection_count = 0
    filament_removed_count = 0
    filament_removed_length = zero(T)
    if criterion(cache) === NoReconnections()
        return (; reconnection_count, filament_removed_count, filament_removed_length,)
    end
    crit = criterion(cache)
    Ls = periods(cache)
    @timeit to "find_reconnection_candidates!" begin
        candidates = find_reconnection_candidates!(cache, fs; to)
    end
    @timeit to "iterate over candidates" for candidate ∈ candidates
        candidate === nothing && continue
        info = should_reconnect(crit, candidate; periods = Ls)
        info === nothing && continue
        info, candidate = find_better_candidates(info, candidate) do other_candidate
            should_reconnect(crit, other_candidate; periods = Ls)
        end
        (; a, b,) = candidate
        @timeit to "reconnect" if a.f === b.f
            # Reconnect filament with itself => split filament into two
            @assert a.i ≠ b.i
            nremoved, Lremoved = reconnect_with_itself!(callback, fs, a.f, a.i, b.i, info)
            filament_removed_count += nremoved
            filament_removed_length += Lremoved
            invalidate_candidates!(cache, a.f)
        else
            # Reconnect two different filaments => merge them into one
            reconnect_with_other!(callback, fs, a.f, b.f, a.i, b.i, info)
            invalidate_candidates!(cache, a.f, b.f)
        end
        reconnection_count += 1
    end
    (; reconnection_count, filament_removed_count, filament_removed_length,)
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
        c′ = ReconnectionCandidate(x′, y′)
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

# Find the index of a filament in a vector of filaments.
# Returns `nothing` if the filament is not found.
# This can happen (but should be very rare) if the filament actually disappeared from `fs`
# after a previous reconnection.
find_filament_index(fs, f) = findfirst(g -> g === f, fs)

# This function reconnects a filament `f` with itself, resulting in two filaments.
# Note that it is possible for one or the two resulting filaments to be too small (in number
# of nodes) to be represented by the discretisation method, in which case they are discarded.
# Hence, this function returns a tuple (nremoved, Lremoved) with:
# - nremoved: number of discarded filaments (0, 1 or 2);
# - Lremoved: total length of discarded filaments, using the straight filament approximation.
function reconnect_with_itself!(callback::F, fs, f, i, j, info) where {F}
    T = number_type(f)
    @assert T <: AbstractFloat

    n = find_filament_index(fs, f)
    n === nothing && return nothing

    gs = split!(f, i, j; p⃗ = info.p⃗)  # = (g₁, g₂)
    @assert length(gs) == 2

    # Determine whether to keep each new filament.
    # We discard a filament if it can't be represented with the chosen discretisation method
    # (typically if the number of nodes is too small)
    keep = map(g -> check_nodes(Bool, g), gs)
    m = findfirst(keep)  # index of first filament to keep (either 1, 2 or nothing)

    if m === nothing
        # We discard the two filaments and remove the original filament from fs.
        popat!(fs, n)
        callback(f, n, :removed)
        nremoved = 2
        Lremoved = sum(filament_length, gs)::T  # note: we use the straight segment approximation to compute lengths (no quadratures)
        return (nremoved, Lremoved)
    end

    # Replace the original filament by the first accepted filament.
    let g = gs[m]
        update_coefficients!(g)
        fs[n] = g
        callback(g, n, :modified)
    end

    if count(keep) == 2
        # If the two new filaments were accepted, append the second one at the end of `fs`.
        @assert m == 1
        let g = gs[2]
            update_coefficients!(g)
            push!(fs, g)
            callback(g, lastindex(fs), :appended)
        end
        nremoved = 0
        Lremoved = zero(T)
    else
        # Otherwise, return the length of the removed filament (which is the first one).
        @assert m == 2
        nremoved = 1
        Lremoved = filament_length(gs[1])::T
    end

    (nremoved, Lremoved)
end

function reconnect_with_other!(callback::F, fs, f, g, i, j, info) where {F}
    @assert f !== g
    nf = find_filament_index(fs, f)
    ng = find_filament_index(fs, g)
    (nf === nothing || ng === nothing) && return nothing  # at least one of the filaments is not present in `fs`

    h = merge!(f, g, i, j; p⃗ = info.p⃗)  # filaments are merged onto `h`
    update_coefficients!(h)

    # Replace `f` by the new filament.
    fs[nf] = h
    callback(h, nf, :modified)

    # Remove `g` from the list of filaments.
    popat!(fs, ng)
    callback(g, ng, :removed)

    nothing
end

reconnect!(cache::AbstractReconnectionCache, args...; kws...) =
    reconnect!(Returns(nothing), cache, args...; kws...)

end
