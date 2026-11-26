# BiotSavart

```@meta
CurrentModule = VortexPasta.BiotSavart
CollapsedDocStrings = false
```

```@docs
BiotSavart
```

## Biot–Savart parameters

### Setting the parameters

```@docs
ParamsBiotSavart
```

### Accessing parameters

```@docs
circulation
periods
domain_is_periodic
```

### Other convenience functions

```@docs
kelvin_wave_period
```

## Ewald splitting kernels

```@docs
AbstractEwaldSplitting
PeriodicSplitting
GaussianSplitting
KaiserBesselSplitting
NoSplitting
```

## Computing Biot–Savart interactions

One can compute velocity (and optionally streamfunction) values on filaments using the
following steps:

1. Define Biot–Savart parameters ([`ParamsBiotSavart`](@ref))
2. Initialise Biot–Savart cache ([`init_cache`](@ref))
3. Compute using existing filaments ([`velocity_on_nodes!`](@ref) or [`compute_on_nodes!`](@ref))

In the last step, filaments are generally passed as a vector `fs` of `n` filaments, `fs = [f₁, f₂, ..., fₙ]`,
where each element `fᵢ` is a [`ClosedFilament`](@ref).

```@docs
init_cache
velocity_on_nodes
velocity_on_nodes!
compute_on_nodes!
```

## Point vortex systems in two dimensions

Biot–Savart interactions in two-dimensional point vortex systems are computed using the same
steps as in 3D with minor differences: 

1. Parameters are defined in [`ParamsPointVortex`](@ref) instead of [`ParamsBiotSavart`](@ref).
2. Instead of a filament vector `fs` one passes a [`PointVortexState`](@ref), which contains
   a vector of vortex points and a vector of signs (`UP` or `DOWN`), to
   [`compute_on_nodes!`](@ref) and [`velocity_on_nodes!`](@ref).
3. The outputs of `compute_on_nodes!` are two-element vectors for `velocity` and scalar
   values for `streamfunction`.

**Example**

```julia
using VortexPasta.BiotSavart
params = ParamsPointVortex(; Γ = 1.0, ...)
cache = BiotSavart.init_cache(params)
points = [Vec2(0.3, 0.1), Vec2(π/3, 2π/3), ...]
signs = [UP, DOWN, ...]  # this is a Vector{SignValue}
vortices = PointVortexState(points, signs)
fields = (
    velocity = similar(vortices.points),
    streamfunction = similar(vortices.points, Float64),  # these are scalar values!
)
compute_on_nodes!(fields, cache, vortices)
```

```@docs
ParamsPointVortex
PointVortexState
SignValue
```

## Short-range interactions

### Backends

```@docs
ShortRangeBackend
NaiveShortRangeBackend
CellListsBackend
```

### Internals

```@docs
ShortRangeCache
max_cutoff_distance
init_cache_short
add_pair_interactions!
foreach_charge
foreach_pair
nearby_charges
```

## Long-range interactions

### Backends

```@docs
LongRangeBackend
NonuniformFFTsBackend
ExactSumBackend
```

### Accessing long-range fields

It may be useful to access computed fields (vorticity, velocity, ...) in Fourier space.
For this, one can use the unexported [`get_longrange_field_fourier`](@ref) function.

```@docs
get_longrange_field_fourier
```

### Internals

```@docs
PointData
LongRangeCache
NullLongRangeCache
LongRangeCacheState
init_cache_long
expected_period
set_num_points!
add_point_charges!
has_real_to_complex
compute_vorticity_fourier!
compute_streamfunction_fourier!
compute_velocity_fourier!
to_coarse_grained_vorticity!
interpolate_to_physical!
transform_to_fourier!
similar(::LongRangeCache, ::Dims{3})
```

## KernelAbstractions utils

[KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) is
used to write generic code which works on CPUs and different kinds of GPUs.

```@docs
ka_generate_kernel
KernelAbstractions.get_backend
KernelAbstractions.device
activate_device!
```

## Internals

```@docs
AbstractBackend
BiotSavartCache
PointVortexCache
HostVector
ZetaHat
background_vorticity_correction!
process_point_charges!
copy_output_values_on_nodes!
```
