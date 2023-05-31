using ..BasicTypes: Zero

"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`BasedOnDistance`](@ref): reconnects filament segments which are closer than a critical distance.
"""
abstract type ReconnectionCriterion end

"""
    distance(crit::ReconnectionCriterion) -> Real

Return the critical distance associated to the reconnection criterion.
"""
function distance end

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = Zero()

"""
    BasedOnDistance <: ReconnectionCriterion
    BasedOnDistance(d_crit; cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

The keyword argument `cos_max` disables the reconnection of nearly parallel segments.
The default value `cos_max = 0.97` disables reconnections when the angle between lines
is ``θ < \\arccos(0.97) ≈ 14°``.
"""
struct BasedOnDistance <: ReconnectionCriterion
    dist    :: Float64
    cos_max :: Float64

    function BasedOnDistance(dist; cos_max = 0.97)
        new(dist, cos_max)
    end
end

distance(c::BasedOnDistance) = c.dist
