# Reconnections

```@meta
CurrentModule = VortexPasta.Reconnections
CollapsedDocStrings = false
```

```@docs
Reconnections
```

Reconnections are generally performed by choosing a [reconnection criterion](@ref Reconnection-criteria) and then calling [`reconnect!`](@ref).

## Functions

```@docs
init_cache
reconnect!
```

## Reconnection criteria

```@docs
ReconnectionCriterion
NoReconnections
ReconnectBasedOnDistance
ReconnectFast
```

## Internals

```@docs
should_reconnect
```
