# BiotSavart

```@meta
CurrentModule = VortexPasta.BiotSavart
CollapsedDocStrings = true
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
process_point_charges!
add_short_range_fields!
local_self_induced_velocity
local_self_induced
nearby_charges
```

## Long-range interactions

### Backends

```@docs
LongRangeBackend
NonuniformFFTsBackend
FINUFFTBackend
CuFINUFFTBackend
ExactSumBackend
```

### Internals

```@docs
LongRangeCache
NullLongRangeCache
init_cache_long
expected_period
folding_limits
set_num_points!
add_point_charges!
has_real_to_complex
compute_vorticity_fourier!
to_smoothed_streamfunction!
to_smoothed_velocity!
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
```

## Internals

```@docs
BiotSavartCache
background_vorticity_correction!
```
