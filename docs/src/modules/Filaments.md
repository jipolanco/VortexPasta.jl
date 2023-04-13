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
discretisation_method
knots
knotlims
Filaments.points
Base.getindex
Base.setindex!
update_coefficients!
normalise_derivatives
normalise_derivatives!
```

## Local discretisations

```@docs
LocalDiscretisationMethod
FiniteDiffMethod
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
CubicSplineMethod
ClosedSplineFilament
```

## Plotting

```@docs
filamentplot!
filamentplot
```

## Index

```@index
```
