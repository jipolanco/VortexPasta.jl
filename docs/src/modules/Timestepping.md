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

### Multirate Runge–Kutta schemes

These schemes are completely explicit, but use different timesteps (and
different RK schemes) for the slow and fast dynamics. They are represented by
an *outer* scheme of order ``n`` with a "large" timestep ``Δt`` for the slowly
evolving terms (which are also expensive to compute), coupled to an *inner*
scheme (typically of order ``n - 1``) with a "small" timestep ``Δt/M``.
The implemented schemes are those described in [Sandu, SIAM J. Numer. Anal. 57 (2019)](https://doi.org/10.1137/18M1205492).

```@docs
MultirateScheme
MultirateMidpoint
SanduMRI33a
SanduMRI45a
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
