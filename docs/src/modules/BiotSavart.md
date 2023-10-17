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
compute_on_nodes!
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
init_cache_short
set_filaments!(c::ShortRangeCache, fs)
add_short_range_velocity!
local_self_induced
local_self_induced_velocity
nearby_segments(c::ShortRangeCache, x⃗)
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
compute_vorticity_fourier!
to_smoothed_streamfunction!
to_smoothed_velocity!
interpolate_to_physical!
transform_to_fourier!
```

## Index

```@index
Pages = ["BiotSavart.md"]
```
