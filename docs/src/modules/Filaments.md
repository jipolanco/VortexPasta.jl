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
minimum_knot_increment
maximum_knot_increment
nodes
Base.getindex
Base.setindex!
update_coefficients!
fold_periodic!
normalise_derivatives
normalise_derivatives!
integrate
find_min_distance
check_nodes
```

## Refinement

```@docs
RefinementCriterion
NoRefinement
BasedOnCurvature
refine!
```

## Reconnections

Reconnections are generally performed by calling [`reconnect!`](@ref) and
according to a chosen reconnection criterion.

```@docs
ReconnectionCriterion
NoReconnections
BasedOnDistance
reconnect!
reconnect_self!
reconnect_other!
split!
merge!
```

## Geometric quantities

The following types are provided as a convenient way of evaluating scalar and
vector quantities of interest along filaments.

```@docs
GeometricQuantity
UnitTangent
CurvatureVector
CurvatureScalar
CurvatureBinormal
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
Pages = ["Filaments.md"]
```
