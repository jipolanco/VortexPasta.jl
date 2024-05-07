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
                   deperiodise_separation,
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
    ) -> Int

Perform filament reconnections according to chosen criterion.

Returns the number of performed reconnections.

Note that, when a filament self-reconnects, this creates new filaments, which
are appended at the end of `fs`.

Moreover, this function will remove reconnected filaments if their number of nodes is too small
(typically ``< 3``, see [`check_nodes`](@ref)).

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
    number_of_reconnections = 0
    criterion(cache) === NoReconnections() && return number_of_reconnections
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
            reconnect_with_itself!(callback, fs, a.f, a.i, b.i, info)
            invalidate_candidates!(cache, a.f)
        else
            # Reconnect two different filaments => merge them into one
            reconnect_with_other!(callback, fs, a.f, b.f, a.i, b.i, info)
            invalidate_candidates!(cache, a.f, b.f)
        end
        number_of_reconnections += 1
    end
    number_of_reconnections
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

function reconnect_with_itself!(callback::F, fs, f, i, j, info) where {F}
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
        return nothing
    end

    # Replace the original filament by the first accepted filament.
    let g = gs[m]
        update_coefficients!(g)
        fs[n] = g
        callback(g, n, :modified)
    end

    # If the two new filaments were accepted, append the second one at the end of `fs`.
    if count(keep) == 2
        @assert m == 1
        let g = gs[2]
            update_coefficients!(g)
            push!(fs, g)
            callback(g, lastindex(fs), :appended)
        end
    end

    nothing
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
