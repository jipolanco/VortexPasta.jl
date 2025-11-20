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
autotune
```

### Accessing parameters

```@docs
circulation
periods
domain_is_periodic
```

## Exported functions

```@docs
init_cache
velocity_on_nodes
velocity_on_nodes!
compute_on_nodes!
```

## Other convenience functions

```@docs
kelvin_wave_period
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
local_self_induced_velocity
local_self_induced
foreach_charge
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
folding_limits
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
```

## Internals

```@docs
AbstractBackend
BiotSavartCache
background_vorticity_correction!
process_point_charges!
copy_output_values_on_nodes!
```
