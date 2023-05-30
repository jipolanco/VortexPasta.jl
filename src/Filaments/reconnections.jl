"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`BasedOnDistance`](@ref): reconnects filament segments which are closer than a critical distance.
"""
abstract type ReconnectionCriterion end

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

"""
    BasedOnDistance <: ReconnectionCriterion
    BasedOnDistance(d_crit)

Reconnects filament segments which are at a distance `d < d_crit`.
"""
struct BasedOnDistance <: ReconnectionCriterion
    dist :: Float64
end
