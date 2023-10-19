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

should_reconnect(c::ReconnectionCriterion, a::Segment, b::Segment; kws...) =
    should_reconnect(c, a.f, b.f, a.i, b.i; kws...)

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = nothing
should_reconnect(::NoReconnections, fx, fy, i, j; kws...) = nothing

"""
    ReconnectBasedOnDistance <: ReconnectionCriterion
    ReconnectBasedOnDistance(d_crit; decrease_length = true, cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

# Optional keyword arguments

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
    decrease_length :: Bool

    function ReconnectBasedOnDistance(dist; cos_max = 0.97, decrease_length = true)
        new(dist, dist^2, cos_max, cos_max^2, decrease_length)
    end
end

distance(c::ReconnectBasedOnDistance) = c.dist

function should_reconnect(
        c::ReconnectBasedOnDistance, fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
    )
    (; dist_sq, cos_max_sq, decrease_length,) = c

    min_dist = find_min_distance(fx, fy, i, j; periods)
    (; d⃗, p⃗, ζx, ζy,) = min_dist
    d² = sum(abs2, d⃗)
    d² > dist_sq && return nothing  # don't reconnect

    # Make sure that reconnections reduce the total length (makes sense energetically for vortices).
    if decrease_length
        length_before = norm(fx[i + 1] - fx[i]) + norm(fy[j + 1] - fy[j])
        length_after = norm(fy[j + 1] - fx[i] - p⃗) + norm(fx[i + 1] - fy[j] + p⃗)
        length_after > length_before && return nothing
    end

    X′ = fx(i, ζx, Derivative(1))
    Y′ = fy(j, ζy, Derivative(1))

    # For now, only return the output of find_min_distance + d² if segments should
    # reconnect.
    success = (; min_dist..., d²,)

    xy = X′ ⋅ Y′
    xy < 0 && return success  # always reconnect antiparallel vortices

    cos² = (xy * xy) / (sum(abs2, X′) * sum(abs2, Y′))
    cos² < cos_max_sq ? success : nothing
end
