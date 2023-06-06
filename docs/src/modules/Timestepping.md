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
TimeInfo
```

## Exported functions

```@docs
init
solve!
step!
```

## Timesteppers

The following timesteppers are exported.
When possible, names are the same as those used by [DifferentialEquations.jl solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).

```@docs
Euler
RK4
SSPRK33
DP5
```

## Adaptivity

```@docs
AdaptivityCriterion
NoAdaptivity
AdaptBasedOnSegmentLength
AdaptBasedOnVelocity
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
