# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

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
