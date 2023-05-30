"""
    ReconnectionCriterion

Abstract type describing a criterion for vortex reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`BasedOnDistance`](@ref): reconnects vortex segments which are closer than a critical distance.
"""
abstract type ReconnectionCriterion end

"""
    NoReconnections <: ReconnectionCriterion

Used to disable vortex reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

"""
    BasedOnDistance <: ReconnectionCriterion
    BasedOnDistance(d_crit)

Reconnects vortex segments which are at a distance `d < d_crit`.
"""
struct BasedOnDistance <: ReconnectionCriterion
    dist :: Float64
end
