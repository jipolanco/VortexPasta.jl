# Timestepping

```@meta
CurrentModule = VortexPasta.Timestepping
```

```@docs
Timestepping
```

## Types

```@docs
VortexFilamentProblem
VortexFilamentSolver
```

## Exported functions

```@docs
init
solve!
step!
```

## Timesteppers

The following timesteppers are exported.
For convenience, names are the same as those used in the [DifferentialEquations.jl](@ref) ecosystem.

```@docs
Euler
RK4
SSPRK33
```

## Internals

```@docs
ExplicitTemporalScheme
TemporalSchemeCache
```

## Index

```@index
Pages = ["Timestepping.md"]
```
