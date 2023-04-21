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

## Functions

```@docs
init_cache
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
init_cache_long
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
