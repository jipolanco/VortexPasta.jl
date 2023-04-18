# BiotSavart

```@meta
CurrentModule = VortexFilamentEwald.BiotSavart
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
CellListMapBackend
```

### Internals

```@docs
ShortRangeCache
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
copy_interpolated_data!
```

## Index

```@index
```
