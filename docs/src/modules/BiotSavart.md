# BiotSavart

```@meta
CurrentModule = VortexPasta.BiotSavart
```

```@docs
BiotSavart
```

## Types

```@docs
ParamsBiotSavart
BiotSavartCache
```

## Exported functions

```@docs
init_cache
velocity_on_nodes!
```

## Short-range interactions

### Backends

```@docs
ShortRangeBackend
NaiveShortRangeBackend
```

### Internals

```@docs
ShortRangeCache
init_cache_short
short_range_velocity
local_self_induced_velocity
add_short_range_velocity_self!
add_short_range_velocity_other!
```

## Long-range interactions

### Backends

```@docs
LongRangeBackend
FINUFFTBackend
ExactSumBackend
```

### Internals

```@docs
LongRangeCache
NullLongRangeCache
init_cache_long
expected_period
folding_limits
reset_fields!
set_num_points!
add_pointcharge!
add_point!
transform_to_fourier!
to_filtered_velocity!
interpolate_to_physical!
long_range_velocity_fourier!
long_range_velocity_physical!
add_long_range_velocity!
```

## Index

```@index
Pages = ["BiotSavart.md"]
```
