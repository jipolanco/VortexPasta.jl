"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`ReconnectBasedOnDistance`](@ref): reconnects filament segments which are closer than a
  critical distance.
"""
abstract type ReconnectionCriterion end

"""
    should_reconnect(
        c::ReconnectionCriterion,
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

should_reconnect(c::ReconnectionCriterion, candidate::ReconnectionCandidate; kws...) =
    should_reconnect(c, candidate.a, candidate.b; kws...)
should_reconnect(c::ReconnectionCriterion, a::Segment, b::Segment; kws...) =
    should_reconnect(c, a.f, b.f, a.i, b.i; kws...)

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = nothing
should_reconnect(::NoReconnections, fx, fy, i, j; kws...) = nothing
max_passes(::NoReconnections) = 0

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
    max_passes  :: Int
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

function Base.show(io::IO, c::ReconnectBasedOnDistance)
    print(io, "ReconnectBasedOnDistance($(c.dist); cos_max = $(c.cos_max), max_passes = $(c.max_passes), use_velocity = $(c.use_velocity), decrease_length = $(c.decrease_length))")
end

function should_reconnect(
        c::ReconnectBasedOnDistance, fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
        can_return_nothing = Val(true),  # this is set to false to obtain the return type when constructing a ReconnectionCache
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
