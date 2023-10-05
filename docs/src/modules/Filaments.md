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
```

## Initialisation

```@docs
Filaments.init
Filaments.from_vector_field
```

## Obtaining information

```@docs
discretisation_method
knots
knotlims
end_to_end_offset
minimum_knot_increment
maximum_knot_increment
nodes
Filaments.distance_to_field
```

## Modifying filaments

```@docs
update_coefficients!
change_offset
fold_periodic!
redistribute_nodes!
```

## Other functions

```@docs
Base.getindex(f::AbstractFilament, i::Int)
Base.setindex!(f::AbstractFilament, v, i::Int)
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
RefineBasedOnSegmentLength
RefineBasedOnCurvature
refine!
insert_node!
remove_node!
update_after_changing_nodes!
```

## Reconnections

Reconnections are generally performed by calling [`reconnect!`](@ref) and
according to a chosen reconnection criterion.

```@docs
ReconnectionCriterion
NoReconnections
ReconnectBasedOnDistance
should_reconnect
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
```

## Segments

### Segment iterators

```@docs
SegmentIterator
segments
length(::SegmentIterator)
eachindex(::SegmentIterator)
```

### Single segments

```@docs
Segment
midpoint
Filaments.segment_length
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
