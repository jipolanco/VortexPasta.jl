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
```

## Reconnection criteria

```@docs
ReconnectionCriterion
NoReconnections
ReconnectBasedOnDistance
```

## Internals

```@docs
should_reconnect
```

## Index

```@index
Pages = ["Reconnections.md"]
```
