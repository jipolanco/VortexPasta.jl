# Filaments

```@meta
CurrentModule = VortexPasta.Filaments
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

## Segments

```@docs
Segments
segments
length(::Segments)
eachindex(::Segments)
```

## Plotting

```@docs
filamentplot!
filamentplot
```

## Index

```@index
```
