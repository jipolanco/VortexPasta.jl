# Diagnostics

```@meta
CurrentModule = VortexPasta.Diagnostics
CollapsedDocStrings = true
```

```@docs
Diagnostics
```

## Kinetic energy

```@docs
kinetic_energy
kinetic_energy_from_streamfunction
kinetic_energy_nonperiodic
```

## Helicity

```@docs
helicity
```

## Filament length

See [`filament_length`](@ref) in the [`Filaments`](@ref) module.
For convenience, `filament_length` is re-exported by `Diagnostics`, meaning that one can do:

```julia
using VortexPasta.Diagnostics
filament_length(...)
```

without needing to import `VortexPasta.Filaments`.

## Stretching rate

```@docs
stretching_rate
```

## Vortex impulse

```@docs
vortex_impulse
```

## Energy spectrum

```@docs
energy_spectrum!
energy_spectrum
init_energy_spectrum
```
