# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

## [0.24.8] - 2024-11-14

### Changed

- Reconnections now use instantaneous filament velocity information by default.
  For now, this is used to discard reconnections between points which are
  getting away from each other (according to their instantaneous velocities).
  This behaviour can be disabled by passing `use_velocity = false` to the
  `ReconnectBasedOnDistance` criterion.

- Remove `BasicTypes` module.
  Definitions have been moved onto the `Filament` module (`Vec3`, `Derivative`)
  and onto the new `Constants` (`Zero`, `Infinity`, `∞`) and `Containers`
  (`VectorOfVectors`) modules.

### Added

- Allow passing `p::ParamsBiotSavart` to `Diagnostics.kinetic_energy*` and `Diagnostics.helicity`.
  This can be more convenient that explicitly passing the circulation `Γ` and maybe the periods `Ls` when needed.

## [0.24.7] - 2024-10-29

### Changed

- Move `FINUFFTBackend` to package extension. This means that one now needs to
  explicitly load FINUFFT.jl to use this backend (#20).

## [0.24.6] - 2024-10-29

### Added

- Add support for AMD GPUs via [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl).
  This requires using the `NonuniformFFTsBackend` for computing long-range
  interactions, and setting the device to `ROCBackend()`.

## [0.24.5] - 2024-09-20

### Added

- Allow using GPU version of NonuniformFFTs.jl for computing long-range interactions.
  For this, we now require NonuniformFFTs.jl v0.5.1.

## [0.24.4] - 2024-09-12

### Added

- Add `VortexPastaThreadPinningExt` extension. When ThreadPinning.jl is loaded
  and `pinthreads` has been used (as recommended when using threads), the
  extension automatically unpins the threads before calls to the FINUFFT
  library (based on OpenMP; seems to conflict with ThreadPinning), and restores
  the pinning afterwards.
  This should improve performance when using the `FINUFFTBackend` with multiple
  threads and ThreadPinning.
  Note that this doesn't affect the GPU-based `CuFINUFFTBackend` or other
  long-range backends.

## [0.24.3] - 2024-09-09

### Added

- Add `BiotSavart.autotune` function, which attempts to determine the optimal
  Ewald parameters based on a given vortex configuration and accuracy parameter $β$.

### Changed

- `CuFINUFFTBackend`: set default upsampling factor to 1.25 (instead of 2.0),
  which seems to be faster for relevant accuracy levels and should use less memory.

- `ParamsBiotSavart`: `rcut` must now be explicitly passed and no longer has a default value.
  (Breaking change, but minor.)

## [0.24.2] - 2024-09-04

### Added

- Performance: long-range computations with GPU backends are now asynchronous,
  and are done while the CPUs compute short-range interactions.

## [0.24.1] - 2024-09-02

### Fixed

- Fix memory "leak" of FINUFFT backends.
  One must explicitly call the `(cu)finufft_destroy!` functions to free the
  memory associated to a plan (FINUFFT.jl could do this for us...).

## [0.24.0] - 2024-08-30

### Added

- Add experimental GPU backend (`CuFINUFFTBackend`) for long-range computations.
  It is based on the CUDA version of the FINUFFT library, meaning that only
  Nvidia GPUs are supported for now.
  It requires FINUFFT 2.3.0-rc1, which is not currently bundled with the
  FINUFFT.jl interface. That is, one must manually compile the FINUFFT
  binaries and link them to the Julia interface.

### Fixed

- Fix `map(f, us::VectorOfVectors)` when the function `f` returns scalar values.
  For instance, this fixes `map(length, us)` when one wants to get the length
  of each contained vector.

## [0.23.1] - 2024-07-09

### Fixed

- Fix computation of squared distance in short-range interactions.
  As a result of this bug, the short-range part was including too many
  interactions, and not correctly taking into account the cut-off distance.
  This means that results were *too precise* compared to the wanted precision.
  This mainly affected the `NaiveShortRangeBackend` (where *all* interactions
  were computed regardless of the given cut-off distance), while the
  `CellListsBackend` was only weakly affected.
  This bug was introduced in v0.19.1.

## [0.23.0] - 2024-07-08

### Changed

- Timestepping: the velocity and streamfunction fields are now represented as
  `ClosedFilament`s. This can be convenient for interpolating their values
  in-between discretisation points.

- Change order of arguments in `Diagnostics.kinetic_energy*` functions, to make
  it more consistent with other diagnostics functions.
  This only applies to the variants which explicitly take a list of filaments and
  streamfunction values.

### Added

- Timestepping: add optional `stretching_velocity` argument to `init`.
  It allows to impose an external velocity on the filaments, parallel to their
  local curvature vector, which works as an artificial stretching (or
  shrinking) mechanism.

## [0.22.2] - 2024-06-24

### Added

- Add `Diagnostics.stretching_rate`, which allows to compute the stretching
  rate of vortex filaments based on their instantaneous velocity.

## [0.22.1] - 2024-06-20

### Added

- Count decrease of vortex length due to reconnection events during a simulation.
  Results are stored in the new `iter.stats.reconnection_length_loss` field.

## [0.22.0] - 2024-06-17

### Added

- Better support for non-zero mean vorticity in the domain.
  We correct the streamfunction so that results do not depend on the splitting parameter $α$.
  Note that the velocity doesn't need to be corrected, as it's not affected by the background vorticity.

## [0.21.3] - 2024-06-04

### Fixed

- Fix wrong length of removed filaments in simulation stats.

## [0.21.2] - 2024-06-03

### Fixed

- Fix bug due to typo when counting removed vortices after reconnections.

## [0.21.1] - 2024-06-03

### Fixed

- Fix possible assertion error in reconnections when one filament must be
  removed after a self-reconnection.

## [0.21.0] - 2024-06-03

### Added

- Add `stats` field to `VortexFilamentSolver`. This allows to know in particular
  the total number of reconnections since the start of a simulation. It also
  stores information about the total number of filaments removed from the
  simulations (because they had not enough nodes) and their total length.

## [0.20.1] - 2024-05-16

### Changed

- `AdaptBasedOnVelocity`: try to reduce number of rejected timesteps.
  For this, we now avoid increasing the *a priori* `dt` too abruptly.

## [0.20.0] - 2024-05-15

### Fixed

- Fix thread safety issue when iterating over elements of a `PeriodicCellList`.
  Basically, when multiple threads iterated over the same `PeriodicCellList`
  (when computing short-range interactions), all threads modified the state of
  the iterator at the same time leading to wrong results.
  We fix this by removing mutable state in `PeriodicCellList` (the old
  `iterators` vector), so that all threads can iterate at the same time.

### Changed

- Timestepping: make sure the `AdaptBasedOnVelocity` criterion is actually satisfied a posteriori.
  If it's not (i.e. if the actual node displacements are too large within
  a timestep), then the displacements are rejected and recomputed with a smaller timestep.
  This can mean longer but more accurate simulations.

## [0.19.4] - 2024-05-13

### Fixed

- Yet another fix for error when writing VTKHDF files with periodic wrapping.
  The `unsafe_trunc` function was incorrectly used, where `round` should be used instead.

## [0.19.3] - 2024-05-13

### Fixed

- Relax assertion condition when writing VTKHDF files with periodic wrapping.
  This avoids possible assertion errors.

## [0.19.2] - 2024-05-13

### Changed

- Accept MakieCore v0.8 (required to use Makie v0.21).

## [0.19.1] - 2024-05-06

### Changed

- Optimise short-range computations.
  One of the main costs of the short-range part is the computation of `erfc`.
  We now use the implementation in the VectorizationBase.jl package.
  We also try to explicitly vectorise (using SIMD) all operations.

## [0.19.0] - 2024-04-24

### Added

- Define `empty!` for `TimeSeriesFile`. Useful when using this within a callback in timestepping.

- In `Filaments.from_vector_field`, make it easier (and document how) to define
  vortex lines from the *curl* of a vector field.
  This can be convenient when one knows the velocity field, and is too lazy to
  analytically derive the corresponding vorticity field.
  The vorticity field is obtained using automatic differentiation.

### Changed

- Parallelise some more serial operations.

- Slightly modify parallelisation of short-range computations.
  Parallelisation is now done at the level of the list of filaments instead of
  at the level of each filament, which makes sense when there are many filaments.

- Change integration limits of `AdaptiveTanhSinh` quadratures.
  This is to ensure that we never evaluate the integrand on endpoints (in case
  there are singularities there).
  Also, we now use different (smaller) integration limits for `Float32`.

- Improve precision of `Filaments.from_vector_field` near the point where the
  filament is closed.

## [0.18.3] - 2024-04-12

### Added

- Add `AdaptiveTanhSinh` quadrature to `Quadratures` module.
  The [tanh-sinh quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature) behaves nicely
  in the presence of endpoint singularities, and can be easily (and efficiently) made adaptive.

- Use `AdaptiveTanhSinh` when `lia_segment_fraction` is enabled in `ParamsBiotSavart`.
  This is to reduce errors when integrating very close to a singularity (however
  at an increased cost, since more function evaluations are generally needed).

## [0.18.2] - 2024-04-11

### Added

- Add `FilamentIO.TimeSeriesFile` for writing JSON time series files which can
  be read by ParaView.

## [0.18.1] - 2024-04-08

### Fixed

- `FilamentIO`: fix possible assertion error when writing VTKHDF files using the
  `periods` argument.

## [0.18.0] - 2024-04-03

### Changed

- Slightly change convention when computing energy spectra:
  * wavenumber bins are now centred at the wavenumber `k` of interest;
  * we now divide the energy spectrum by the wavenumber increment `Δk`, which
    makes sense if we consider the total energy as the integral of the energy
    spectrum. This makes no difference when the domain size is 2π (in which case
    `Δk = 1`).

- Change default quadrature rule in `ParamsBiotSavart` from `GaussLegendre(2)`
  to `GaussLegendre(3)`.

- `FINUFFTBackend`: perform NUFFTs of all vector components at once. This
  improves performance of this backend. For the same precision, performance is
  now comparable to that of `NonuniformFFTsBackend` provided that the FINUFFT
  libraries are compiled from source (see [here](https://github.com/ludvigak/FINUFFT.jl?tab=readme-ov-file#advanced-installation-and-locally-compiling-binaries) for details).

### Added

- Add `longrange_truncate_spherical` option to `ParamsBiotSavart`. This should
  only be used for testing.

- Add optional `callback_vorticity` argument to `BiotSavart.compute_on_nodes!`.
  It gives access to the truncated Fourier coefficients of the vorticity field,
  it is converted onto the coefficients of the long-range velocity or streamfunction.

- Add optional `periods` argument to `filamentplot!`, which allows to "break"
  a curve onto multiple sub-curves which fall inside the periodic domain.

- Similarly, add optional `periods` argument to `FilamentIO.write_vtkhdf`.
  VTKHDF files with broken lines can be read back onto the original "unbroken"
  filaments.

### Removed

- Remove old `regularise_binormal` option, which is almost always a bad idea
  when estimating the locally-induced velocity.

## [0.17.1] - 2024-02-23

### Added

- Allow computing energy spectrum from unsmoothed vorticity in Fourier space.
  Also, make it easier to compute energy spectra with a different resolution (`kmax`) in
  wavenumber space than the resolution used for long-range computations.

### Changed

- Change some internals in the computation of the short-range Biot–Savart
  component to avoid potential precision issues (just in case; not sure these
  issues were really significant).

## [0.17.0] - 2024-02-14

### Changed

- Reduce default precision of NUFFT backends (`FINUFFTBackend`, `NonuniformFFTsBackend`).
  The relative tolerance is reduced from `1e-8` to `1e-6`, which is still
  largely enough in applications (especially since this tolerance doesn't
  account for the Gaussian smoothing proper of Ewald summation, which should further reduce the error).

- Strang splitting: change default scheme used for the "fast" dynamics to `RK4` (used to be `Midpoint`).

- Parallelise LIA-only computations. This only makes a difference either when (1) doing LIA-only simulations,
  or (2) using splitting/multirate/IMEX timestepping schemes which consider the LIA (local) term as the "fast" term.

## [0.16.0] - 2024-02-06

### Added

- Timestepping: add Strang splitting scheme.

- Timestepping: add `LIA` option to `init`, enabling simulations of the local
  induction approximation.

- Timestepping: add `MaximumTimestep` adaptivity criterion. This criterion is
  implicitly applied when other criteria are passed which allow "infinite"
  timesteps (the case of `AdaptBasedOnVelocity`).

- Write filament parametrisation (knots) to VTKHDF files in `write_vtkhdf`.
  This can be disabled by passing `parametrisation = false`.
  Also, if the parametrisation is present in a VTKHDF file, then `read_vtkhdf`
  reads the existent parametrisation instead of computing a new one.
  This can be useful when using non-standard parametrisations.

- Add custom definitions of `view`, `==` and `isapprox` for `PaddedArray`s.

## [0.15.1] - 2024-01-30

### Added

- Add `Diagnostics.helicity` to compute the helicity of a vortex filament configuration.

### Fixed

- Make sure `integrate` returns `Float32` values when working with `Float32`.
  This fixes some type instabilities when using single precision.

- The behaviour of reconnections when periodicity is disabled is now more
  consistent with the behaviour with periodicity enabled. In particular,
  reconnections with small critical distances should work better than before.

## [0.15.0] - 2024-01-26

### Added

- VTKHDF: add convenient syntax for writing geometric quantities (such as `CurvatureVector`).

- Add predefined curve parametrisations.

### Changed

- VTKHDF: more consistent and accurate interpolation of data on filaments when
  refinement is enabled. Instead of using linear interpolation, we use the same
  interpolation method used to represent the filaments themselves.

- Slightly modify energy computation from streamfunction: when using
  quadratures, we now interpolate the tangent component of the streamfunction on nodes,
  instead of interpolating the streamfunction vector and *then* taking the tangent component.

## [0.14.0] - 2024-01-22

### Changed

- VTKHDF files: slightly change internal representation of datasets (use
  `PolyData` instead of `UnstructuredGrid`).
  This doesn't change anything to visualisation, other than ParaView 5.12 is
  now needed.

### Added

- Add `fold_periodic` option to `Timestepping.init`.

## [0.13.1] - 2024-01-22

### Fixed

- Improve reconnection criteria when setting a critical distance which is much
  smaller than the typical segment length.

## [0.13.0] - 2024-01-19

### Changed

- Improve performance of reconnections in periodic domains (using periodic cell lists).

- Optimise periodic padding of multidimensional `PaddedArray`s.

### Fixed

- Changed implementation of periodic cell lists to fix possible memory leak
  when using a large amount of cells (at least $∼512³$), which can occur when
  enabling reconnections with a small reconnection distance in a periodic domain.
  The issue seems to be related to creating a very large array of arrays, which
  is avoided in the new implementation.

## [0.12.1] - 2024-01-18

### Fixed

- `Diagnostics.kinetic_energy` always uses the streamfunction. Also removed
  warnings when calling `kinetic_energy_from_streamfunctions` in non-periodic
  domains, where this definition is totally valid and should be preferred.

- Fix possible memory leak when using `sizehint!` when emptying a `PeriodicCellList`.
  This was mostly visible when using cell lists for reconnections (in which case each cell is usually very small).

## [0.12.0] - 2024-01-12

### Added

- Add alternative energy estimation method adapted for non-periodic domains.

## [0.11.2] - 2024-01-04

### Added

- Add example: unforced Kelvin waves.

### Fixed

- Fix wrong results when using splitting timestepping schemes (IMEX or
  multirate) and setting an `external_velocity`.

## [0.11.1] - 2024-01-03

### Changed

- `NonuniformFFTsBackend` now uses NonuniformFFTs.jl v0.3, which improves performance.

## [0.11.0] - 2023-12-20

### Added

- Add `NonuniformFFTsBackend` for long-range computations.

### Changed

- `NonuniformFFTsBackend` is now the default backend for long-range computations.

## [0.10.0] - 2023-11-27

### Added

- Add optional `affect!` keyword argument to `Timestepping.init`.
  It works similarly to `callback`, but allows modifying filament definitions
  before computing their induced velocities.

## [0.9.0] - 2023-11-24

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

- New [`FourierMethod`](https://jipolanco.github.io/VortexPasta.jl/dev/modules/Filaments/#VortexPasta.Filaments.FourierMethod) filament discretisation method
  for describing filaments using Fourier series (!53).
  Enables accurate estimation of curve derivatives via FFTs, but only for simple settings where discretisation points stay more or less equispaced.
