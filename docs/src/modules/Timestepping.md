# Timestepping

```@meta
CurrentModule = VortexPasta.Timestepping
CollapsedDocStrings = true
```

```@docs
Timestepping
```

## Setting-up a simulation

The usual way of setting-up a simulation is to first create
a [`VortexFilamentProblem`](@ref) and then to call [`init`](@ref) to initialise
a [`VortexFilamentSolver`](@ref):

```julia
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping

fs = [Filaments.init(...) for n ∈ 1:10]          # initialise a set of filaments
params = ParamsBiotSavart(...)                   # set Biot-Savart parameters
tspan = (tmin, tmax)                             # set start and end time of simulation
prob = VortexFilamentProblem(fs, tspan, params)  # create problem
iter = init(prob, RK4(); dt = 0.01, ...)         # initialise VortexFilamentSolver
```

---

```@docs
VortexFilamentProblem
init
VortexFilamentSolver
```

## Running a simulation

There are basically two ways of running a simulation:

1. either by calling [`solve!(iter)`](@ref), which will run the full simulation
   up to the final time `tmax`;

2. or by repeatedly calling [`step!(iter)`](@ref) (for example inside a `for` loop)
   to advance the simulation one timestep at a time.


The second option can seem to be more convenient as it allows to do things like
running analyses or saving snapshots at intermediate stages of the simulation.
However, those things are also easy to do with the first option, by passing
a `callback` to the [`init`](@ref) function.
See [this section](@ref tutorial-vortex-ring-simulation-state) of the vortex
ring tutorial for examples.


```@docs
solve!
step!
```

## Temporal schemes

The following timesteppers are available.
When possible, names are the same as those used by [DifferentialEquations.jl solvers](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/).

```@docs
TemporalScheme
```

### Explicit Runge–Kutta schemes

```@docs
ExplicitScheme
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
ImplicitExplicitScheme
IMEXEuler
Ascher343
KenCarp3
KenCarp4
```

### Splitting schemes

```@docs
SplittingScheme
Strang
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

### Implicit schemes

These implicit schemes are mainly meant to be used as inner schemes when using
multirate methods.

```@docs
ImplicitScheme
CrankNicolson
```

### Determining the fast term

The splitting between fast and slow terms in IMEX and multirate schemes can be done in two different ways.
One may either choose to identify the fast term with the local (LIA) term in Biot–Savart, or with the short-range component of Ewald summation.

```@docs
LocalTerm
ShortRangeTerm
```

## [Adaptivity criteria](@id Adaptivity)

A detailed below, a few temporal adaptivity criteria are available which can be
used as the `adaptivity` argument of [`init`](@ref).

```@docs
AdaptivityCriterion
NoAdaptivity
AdaptBasedOnSegmentLength
AdaptBasedOnVelocity
MaximumTimestep
```

## Internals

```@docs
TemporalSchemeCache
TimeInfo
```
