# Filaments

```@meta
CurrentModule = VortexPasta.Filaments
CollapsedDocStrings = false
```

```@docs
Filaments
```

## Types

```@docs
AbstractFilament
ClosedFilament
Vec3
Derivative
```

## Initialisation

```@docs
Filaments.init
Filaments.from_vector_field
```

## Curve representation

### Discretisation methods

```@docs
DiscretisationMethod
SplineMethod
CubicSplineMethod
QuinticSplineMethod
FiniteDiffMethod
FourierMethod
discretisation_method
```

### Interpolation

```@docs
interpolation_method
Filaments.required_derivatives
Filaments.init_coefficients
Filaments.compute_coefficients!
Filaments.evaluate
LocalInterpolationMethod
HermiteInterpolation
interpolate
```

## Obtaining information

```@docs
knots
knotlims
end_to_end_offset
minimum
maximum
minimum_node_distance
minimum_knot_increment
maximum_knot_increment
nodes
filament_length
Filaments.distance_to_field
Filaments.minimum_nodes
```

## Modifying filaments

```@docs
update_coefficients!
close_filament!
change_offset
fold_periodic!
redistribute_nodes!
split!
merge!
```

## Iterating over filament nodes

```@docs
FilamentChunkIterator
```

## Curve parametrisations

```@docs
AbstractParametrisation
ChordalParametrisation
CentripetalParametrisation
FourierParametrisation
CustomParametrisation
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
number_type
Filaments.curl
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

## Geometric quantities

The following types are provided as a convenient way of evaluating scalar and
vector quantities of interest along filaments.

```@docs
GeometricQuantity
UnitTangent
CurvatureVector
CurvatureScalar
CurvatureBinormal
TorsionScalar
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
