"""
    Reconnections

Module for dealing with the reconnection of filaments.
"""
module Reconnections

export
    ReconnectionCriterion,
    NoReconnections,
    ReconnectBasedOnDistance,
    ReconnectFast,
    reconnect!

using ..Constants: Infinity

using ..Filaments:
    Filaments,
    AbstractFilament,
    ClosedFilament,
    Segment,
    Derivative,
    Vec3,
    segments,
    split!, merge!,
    update_coefficients!,
    find_min_distance,
    filament_length,
    deperiodise_separation,
    number_type,
    check_nodes

using ..CellLists: CellLists

using OhMyThreads: OhMyThreads, tforeach, DynamicScheduler
using TimerOutputs
using LinearAlgebra: ⋅, norm

"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`ReconnectBasedOnDistance`](@ref): reconnects filament segments which are closer than a
  critical distance;

- [`ReconnectFast`](@ref): considers segments as straight lines, enabling to perform all
  required reconnections in a single pass.
"""
abstract type ReconnectionCriterion end

abstract type AbstractReconnectionCache end
distance(c::AbstractReconnectionCache) = distance(criterion(c))
max_passes(c::AbstractReconnectionCache) = max_passes(criterion(c))

"""
    Reconnections.init_cache(
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament},
        Ls::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) -> AbstractReconnectionCache

Initialise reconnection cache.

Required arguments are a reconnection criterion `crit` and the domain dimensions (or *periods*) `Ls`.
The type of cache will vary depending on the chosen reconnection criterion.
"""
function init_cache(
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament},
        Ls::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )
    _init_cache(crit, fs, Ls)
end

struct ReconnectionCandidate{S <: Segment}
    a :: S
    b :: S
    filament_idx_a :: Int  # filament id where segment `a` is located
    filament_idx_b :: Int  # filament id where segment `b` is located
end

include("criteria/no_reconnections.jl")
include("criteria/based_on_distance.jl")
include("criteria/reconnect_fast.jl")

"""
    reconnect!(
        [callback::Function],
        cache::AbstractReconnectionCache,
        fs::AbstractVector{<:AbstractFilament},
        [vs::AbstractVector{<:AbstractVector}],
    ) -> NamedTuple

Perform filament reconnections according to chosen criterion.

Note that, when a filament self-reconnects, this creates new filaments, which
are appended at the end of `fs`.

Moreover, this function will remove reconnected filaments if their number of nodes is too small
(typically ``< 3``, see [`check_nodes`](@ref)).

Optionally, `vs` can be a vector containing all instantaneous filament velocities.
In this case, this will be used to discard reconnections between filaments which are moving
in opposite directions. This is required if the [`ReconnectBasedOnDistance`](@ref) criterion
is used with its option `use_velocity = true` (default if `max_passes = 1`).

## Returns

This function returns a `NamedTuple` with fields:

- `reconnection_count`: **number of reconnections**;

- `reconnection_length_loss`: **decrease of filament length due to reconnections**. This is
  simply an estimate for the difference in filament length before and after a reconnection.
  It is unrelated to (and does not include) filament removal (see below).
  Uses the straight segment approximation;

- `filaments_removed_count`: **number of filaments removed after reconnection**. This happens
  when a filament reconnects with itself, splitting into two new filaments, and one or both
  of these filaments have an insufficient number of discretisation nodes;

- `filaments_removed_length`: **length of removed filaments**. We use a straight segment
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
        fs::AbstractVector{<:AbstractFilament},
        vs::Union{Nothing, AbstractVector{<:AbstractFilament}} = nothing;
        to::TimerOutput = TimerOutput(),
    ) where {F <: Function}
    @timeit to "Reconnect" begin
        T = number_type(fs)
        ret = (;
            reconnection_count = 0,
            reconnection_length_loss = zero(T),  # decrease of vortex length due to reconnections (not due to removals)
            filaments_removed_count = 0,
            filaments_removed_length = zero(T),
            npasses = 0,
        )
        R = typeof(ret)
        nmax = max_passes(cache)
        npasses = 0
        while npasses < nmax
            npasses += 1
            nrec = ret.reconnection_count
            ret = _reconnect_pass!(callback, ret, cache, fs, vs; to)::R
            nrec == ret.reconnection_count && break  # no new reconnections were performed
        end
    end
    ret
end

reconnect!(cache::AbstractReconnectionCache, args...; kws...) =
    reconnect!(Returns(nothing), cache, args...; kws...)

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
    n === nothing && error("filament `f` not found in vector of filaments `fs`")

    gs = split!(f, i, j; p⃗ = info.p⃗)  # = (g₁, g₂)
    @assert length(gs) == 2

    # Determine whether to keep each new filament.
    # We discard a filament if it can't be represented with the chosen discretisation method
    # (typically if the number of nodes is too small)
    keep = map(g -> check_nodes(Bool, g), gs)
    m = findfirst(keep)  # index of first filament to keep (either 1, 2 or nothing)
    for g ∈ gs
        Filaments.close_filament!(g)  # needed before calling filament_length
    end

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
        # Otherwise, return the length of the removed filament.
        @assert count(keep) == 1
        m_removed = 3 - m  # index of removed filament ∈ (1, 2)
        nremoved = 1
        Lremoved = filament_length(gs[m_removed])::T
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

end
