"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = nothing
should_reconnect(::NoReconnections, fx, fy, i, j; kws...) = nothing
max_passes(::NoReconnections) = 0

# Used when reconnections are disabled.
struct NullReconnectionCache <: AbstractReconnectionCache end
_init_cache(::NoReconnections, args...) = NullReconnectionCache()
criterion(::NullReconnectionCache) = NoReconnections()

# Don't do anything; just return `ret_base` for consistency with other criteria.
function _reconnect_pass!(callback::F, ret_base, cache::NullReconnectionCache, args...; kwargs...) where {F}
    ret_base
end
