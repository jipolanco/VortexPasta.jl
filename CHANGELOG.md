# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Changed

- Move reconnections to separate `Reconnections` module.

- Move most of the short-range backend implementation to new `FindNearbySegments` module.

## [0.4.0] - 2023-10-16

### Added

- New [`FourierMethod`](https://jipolanco.pages.in2p3.fr/VortexPasta.jl/modules/Filaments/#VortexPasta.Filaments.FourierMethod) filament discretisation method
  for describing filaments using Fourier series (!53).
  Enables accurate estimation of curve derivatives via FFTs, but only for simple settings where discretisation points stay more or less equispaced.
