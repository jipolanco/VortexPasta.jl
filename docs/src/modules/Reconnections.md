# Reconnections

```@meta
CurrentModule = VortexPasta.Reconnections
```

```@docs
Reconnections
```

Reconnections are generally performed by choosing a [reconnection criterion](@ref Reconnection-criteria) and then calling [`reconnect!`](@ref).

## Functions

```@docs
init_cache
reconnect!
reconnect_self!
reconnect_other!
should_reconnect
```

## Reconnection criteria

```@docs
ReconnectionCriterion
NoReconnections
ReconnectBasedOnDistance
```

## Index

```@index
Pages = ["Reconnections.md"]
```
