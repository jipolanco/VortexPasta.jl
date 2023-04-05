# Filaments

```@meta
CurrentModule = VortexFilamentEwald.Filaments
```

```@docs
Filaments
```

## Types

```@docs
DiscretisationMethod
AbstractFilament
ClosedFilament
PaddedVector
Derivative
Vec3
```

## Functions

```@docs
Filaments.init
nodes
estimate_derivatives!
normalise_derivatives
normalise_derivatives!
derivatives
derivative
```

## Local discretisations

```@docs
LocalDiscretisationMethod
FiniteDiff
ClosedLocalFilament
```

### Local interpolations

```@docs
LocalInterpolationMethod
HermiteInterpolation
interpolate
```

## Global discretisations

```@docs
GlobalDiscretisationMethod
ClosedSplineFilament
```
