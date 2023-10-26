# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Add [`QuinticSplineMethod`] for accurate filament discretisation.

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
