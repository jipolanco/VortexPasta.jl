# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## 0.12.1 - 2024-01-18

### Fixed

- `Diagnostics.kinetic_energy` always uses the streamfunction. Also removed
  warnings when calling `kinetic_energy_from_streamfunctions` in non-periodic
  domains, where this definition is totally valid and should be preferred.

- Fix possible memory leak when using `sizehint!` when emptying a `PeriodicCellList`.
  This was mostly visible when using cell lists for reconnections (in which case each cell is usually very small).

## 0.12.0 - 2024-01-12

### Added

- Add alternative energy estimation method adapted for non-periodic domains.

## 0.11.2 - 2024-01-04

### Added

- Add example: unforced Kelvin waves.

### Fixed

- Fix wrong results when using splitting timestepping schemes (IMEX or
  multirate) and setting an `external_velocity`.

## 0.11.1 - 2024-01-03

### Changed

- `NonuniformFFTsBackend` now uses NonuniformFFTs.jl v0.3, which improves performance.

## 0.11.0 - 2023-12-20

### Added

- Add `NonuniformFFTsBackend` for long-range computations.

### Changed

- `NonuniformFFTsBackend` is now the default backend for long-range computations.

## 0.10.0 - 2023-11-27

### Added

- Add optional `affect!` keyword argument to `Timestepping.init`.
  It works similarly to `callback`, but allows modifying filament definitions
  before computing their induced velocities.

## 0.9.0 - 2023-11-24

### Added

- Add an option to force the vortex filaments by adding an external velocity.
  This is done via the `external_velocity` keyword argument of `Timestepping.init`.

## [0.8.2] - 2023-11-22

### Added

- Parallelise short-range computations using threads.
  Note that the most expensive part of long-range computations (via FINUFFT) are also parallel.
  To use multithreading, start Julia with the `-t N` flag, where `N` is the number of threads.
  See the [docs](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads) for more details.

## [0.8.1] - 2023-11-10

### Changed

- Relax number of allowed nodes for `QuinticSplineMethod` from 7 to 5.

- Try to do all computations with the same float type. This should improve consistency when using precisions other than `Float64` (for instance `Float32`).

## [0.8.0] - 2023-11-08

### Added

- Add optional `lia_segment_limits` parameter, which allows to compute the LIA term over a fraction of the local segments. This can improve accuracy, especially when the filament discretisation is relatively coarse.

- Add optional `fast_term` parameter in `Timestepping`. When using split schemes (such as IMEX or multirate Runge–Kutta), this allows to indicate what is considered to be the "fast" term. More precisely, the fast term can either be identified to the local (LIA) term in Biot–Savart, or to the short-range term in Ewald summation.

## [0.7.1] - 2023-11-06

### Fixed

- Fix knot insertion issue for quintic splines.

## [0.7.0] - 2023-11-03

### Changed

- `ParamsBiotSavart` now takes a single `quadrature` argument instead of separate `quadrature_short` and `quadrature_long` arguments.
  This allows to reuse interpolation operations between short-range and long-range computations, where positions and derivatives are evaluated in-between discretisation points.
  In particular, short-range computations should be faster now.

### Fixed

- `QuinticSplineMethod` now requires filaments with at least 7 nodes, due to how the associated linear system is assembled.

## [0.6.4] - 2023-10-31

### Fixed

- Make sure the number of derivatives of a filament is preserved when calling `similar` and `copy`.

## [0.6.3] - 2023-10-31

### Added

- Enable computation of curve torsion (`TorsionScalar`).
  This is only possible with high-order filament discretisations (`QuinticSplineMethod`, `FourierMethod`).

## [0.6.2] - 2023-10-30

### Changed

- Optimise internal implementation of `QuinticSplineMethod` for better performance.

## [0.6.1] - 2023-10-27

### Changed

- Change definition of `AdaptBasedOnSegmentLength` based on a more physical criterion, namely the Kelvin wave period at a certain scale.

## [0.6.0] - 2023-10-26

### Added

- Add `QuinticSplineMethod` for accurate filament discretisation.

## [0.5.1] - 2023-10-24

### Added

- Add and export `Filaments.minimum_node_distance` function.

## [0.5.0] - 2023-10-23

### Added

- Add a few options to `Filaments.from_vector_field`.

### Changed

- Change default estimation of binormal vector. Previously it was regularised, now it's more physically correct (!58).

- Move reconnections to separate `Reconnections` module (!54).

- Move most of the short-range backend implementation to new `FindNearbySegments` module (!55).

### Fixed

- Fix periodicity effects on reconnections in `Timestepping` module (!57).

- Fix node insertion and curve refinement when using `FiniteDiffMethod`.

## [0.4.0] - 2023-10-16

### Added

- New [`FourierMethod`](https://jipolanco.pages.in2p3.fr/VortexPasta.jl/modules/Filaments/#VortexPasta.Filaments.FourierMethod) filament discretisation method
  for describing filaments using Fourier series (!53).
  Enables accurate estimation of curve derivatives via FFTs, but only for simple settings where discretisation points stay more or less equispaced.
