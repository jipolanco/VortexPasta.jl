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

## Temporal schemes

The following timesteppers are exported.
When possible, names are the same as those used by [DifferentialEquations.jl solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).

### Explicit Runge–Kutta schemes

```@docs
Euler
Midpoint
RK4
SSPRK33
DP5
```

### Implicit-explicit Runge–Kutta (IMEX-RK) schemes

The following schemes treat local interactions implicitly and non-local interactions explicitly.
This should hopefully allow for larger timesteps than fully explicit schemes.

```@docs
IMEXEuler
Ascher343
KenCarp3
KenCarp4
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
ExplicitScheme
TemporalSchemeCache
```

## Index

```@index
Pages = ["Timestepping.md"]
```
