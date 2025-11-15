# Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## Unreleased

### Fixed

- Fix possible issue with `ReconnectFast` when running in serial mode (`nthreads = 1`).

## [0.30.2] - 2025-11-14

### Added

- Add `CellLists.foreach_pair` and `CellLists.foreach_source` as alternative
  iteration methods over pairs. These are faster than the iterator-based approach
  and more amenable to GPU computing (especially `foreach_pair`). These are now
  used in reconnections and some Biot–Savart computations (excluding the SIMD-accelerated
  version of short-range interactions, where things are more complicated).

### Changed

- Simplify `CellLists` implementation. Now cell lists use less memory and are possibly faster
  than before.

## [0.30.1] - 2025-11-07

### Added

- Allow NonuniformFFTs.jl v0.9.

## [0.30.0] - 2025-11-06

### Changed

- `ParamsBiotSavart`: change default value of `quadrature_near_singularity`. Now it is equal
  to `quadrature` (by default `GaussLegendre(3)`), which is much cheaper and doesn't seem to
  considerably change accuracy.

- BiotSavart: change default computation of short-range integrals and the correction to the
  "self-interaction", to avoid computation of `erf` (and also `exp`) when correcting long-range computations.
  The new solution may lead to (slight) catastrophic cancellation, which seems
  to be really negligible when working in double precision (`Float64`), and is only
  slightly visible in single precision (`Float32`), as only one test needed to be updated.
  One can pass `avoid_explicit_erf = false` to `ParamsBiotSavart` to fall back to the old behaviour.

- `ParamsBiotSavart`: remove deprecated `quadrature_short` and `quadrature_long` parameters.

## [0.29.17] - 2025-10-23

### Added

- Add new `filament_nderivs` option to `Timestepping.init`.

### Changed

- Minor improvements to speed of implicit and IMEX timestepping schemes.

- The old `alias_u0` option of `Timestepping.init` is now deprecated and ignored. Filaments
  are now always copied, which was already the default previously.

- Use DocumenterVitepress for documentation.

## [0.29.16] - 2025-09-24

### Fixed

- Timestepping: fix combination of `adaptivity = AdaptBasedOnVelocity(...)` and `load_checkpoint(...; read_time = false)`.

## [0.29.15] - 2025-09-24

### Added

- `Timestepping.load_checkpoint`: add `read_time = false` option, allowing to restart a
  simulation from time `t = 0`.

## [0.29.14] - 2025-09-17

### Added

- Add `BiotSavart.velocity_on_nodes` (without the `!`), which doesn't require a preallocated
  vector of velocities.

## [0.29.13] - 2025-09-15

### Added

- Add `Diagnostics.energy_transfer_matrix`.

- Add `modify_length` option to `FourierBandForcingBS`. If set to `false`, the
  forcing velocity will be adapted so that it doesn't locally modify vortex lengths.

- `Diagnostics.energy_flux`: compute additional `fluxes.vinf` field for estimating non-local
  transfer to very small scales.

- `Diagnostics.energy_flux`: allow passing a vector of wavenumbers `ks` instead
  of the wanted number of wavenumbers `Nk`.

## [0.29.12] - 2025-09-09

### Added

- Timestepping: include `iter.vf` field when using `NormalFluidForcing` and `FourierBandForcing`
  to match the behaviour of `FourierBandForcingBS`.

## [0.29.11] - 2025-09-05

### Fixed

- **Fix infamous "restart" bug.**
  In fact this wasn't really an issue with restarts, but an issue with filament buffers used
  in timestepping (for RK substeps), i.e. the `fc` field of `TemporalSchemeCache`. In fact,
  we never updated the end-to-end offset of these filament buffers, even when the offset of
  the actual filaments changed (typically after reconnections).

## [0.29.10] - 2025-08-29

### Added

- Add `Diagnostics.energy_flux`.

### Changed

- Require Julia 1.11.

## [0.29.9] - 2025-08-28

### Changed

- Reconnections: use the new (parallelised) `CellLists.set_elements!` when choosing the
  `ReconnectFast` criterion.

## [0.29.8] - 2025-08-28

### Changed

- Improve CPU parallelisation of short-range interactions.
  This should improve performance in general cases, especially when the number of nodes per
  filament varies greatly from one filament to another (as can be the case in turbulence
  simulations).

## [0.29.7] - 2025-08-27

### Added

- Parallelise cell list initialisation. This is now used to speed-up short-range Biot--Savart computations.

## [0.29.6] - 2025-08-27

### Added

- `FourierBandForcingBS` now accepts an optional `α′` parameter for including
  a "drift" term which doesn't inject energy at the targeted wavenumbers.

## [0.29.5] - 2025-08-22

### Added

- `ParamsBiotSavart`: add `use_simd` option allowing to disable explicit SIMD implementation
  of short-range computations.

### Changed

- Avoid using VectorizationBase.jl for SIMD, as this package is no longer maintained.
  This package seems to cause compilation issues in certain computing clusters
  (this is the case on an AMD MI300A cluster).
  SIMD is used for computing `erf(x)` in the short-range component of Ewald summation.
  We now use the SIMD.jl package instead.

## [0.29.4] - 2025-07-23

### Added

- Support loading VTKHDF file written with a different discretisation method.
  Previously this failed when there were filaments with too few nodes compared to the
  requirements of the chosen method. For example, if the file was written with
  `CubicSplineMethod`, it may include filaments with 3 or 4 nodes, which is not allowed by
  `QuinticSplineMethod`. In this case, filaments are now removed and a warning is shown.

## [0.29.3] - 2025-07-23

### Fixed

- Fix construction of `SimulationStats` (issue introduced in v0.29.2).

## [0.29.2] - 2025-07-23

### Added

- Timestepping: include total number of reconnection passes in `iter.stats`.

## [0.29.1] - 2025-07-23

### Added

- `ReconnectFast`: support `max_passes` option, which may be particularly useful when running with multiple threads.

## [0.29.0] - 2025-07-22

### Added

- Add alternative `ReconnectFast` reconnection criterion.
  Compared to `ReconnectBasedOnDistance`, this criterion is able to perform all required reconnections at once, instead of needing multiple passes.
  Geometrically, this criterion is less accurate than `ReconnectBasedOnDistance` because it considers filament segments as straight.

- Document and extend `Filaments.minimum_nodes`.

## [0.28.8] - 2025-07-04

### Fixed

- Fix simulations with `reconnect = NoReconnections()` (issue introduced in v0.28.7).

## [0.28.7] - 2025-07-04

### Added

- `ReconnectBasedOnDistance`: add `max_passes` option for performing multiple reconnection passes.

## [0.28.6] - 2025-07-04

### Performance

- Various major performance improvements to reconnections.

### Fixed

- Fix `FourierBandForcingBS` on AMDGPU.

## [0.28.5] - 2025-06-30

### Added

- Add `DissipationBS` dissipation method.
  This is a cheaper alternative to `SmallScaleDissipationBS` which acts at all scales.

## [0.28.4] - 2025-06-27

### Added

- Add `Diagnostics.integral_lengthscale`.

## [0.28.3] - 2025-06-23

### Fixed

- Fix evaluation of physical-space quantities (long-range velocity, ...) when using GPUs.
  Previously, callbacks passed to NonuniformFFTs failed to compile on CUDA.
  This should also fix `SmallScaleDissipationBS` on GPUs.
  (This fixes an issue which appeared in v0.28.0.)

- Fix computation of spectra on GPUs, which failed with a compilation error on
  recent versions (v0.28.0 probably?).

## [0.28.2] - 2025-06-23

### Changed

- Require Makie.jl 0.24 for plotting.

## [0.28.1] - 2025-06-20

### Added

- Add `SmallScaleDissipationBS`, which allows to dissipate energy at small scales
  (**experimental**).

- When using `FourierBandForcingBS`, the `iter.vf` field is now available in
  `VortexFilamentSolver` which includes the forcing term.

## [0.28.0] - 2025-06-18

### Changed

- **BREAKING**: velocity and streamfunction fields contained in `LongRangeCache` are **no
  longer Gaussian-filtered** to avoid loss of information when doing operations non related
  to Ewald's method (such as energy spectra, some forcing methods, ...). Instead, filtering
  is done "on the fly" when interpolating long-range fields to vortex locations. This change
  requires modifications (usually simplifications) to all code that uses these Fourier-space
  fields.

- **BREAKING**: as a result of the item above, the `unfilter` keyword argument of
  `Diagnostics.energy_spectrum` and `Diagnostics.helicity_spectrum` is no longer accepted,
  and the returned results are always "unfiltered" (which is what we usually want).

- (Slightly) **BREAKING**: Vorticity in Fourier space now has the right dimensions
  (previously, it was missing a `Γ/V` factor, which was included later when computing the
  streamfunction or velocity).

### Removed

- Remove `FINUFFTBackend` and `CuFINUFFTBackend`. This is in part because we're
  starting to use the new callback feature of NonuniformFFTs.jl to avoid extra
  memory copies when computing NUFFTs. This will be useful to store the
  non-smoothed streamfunction andvelocity fields in Fourier space, while
  transfoming the smoothed (long-range) fields to physical space.

## [0.27.10] - 2025-06-16

### Added

- Allow NonuniformFFTs.jl 0.8.

## [0.27.9] - 2025-06-16

### Changed

- The internal `BiotSavart.get_parameters(::LongRangeCache)` function now
  returns all Biot–Savart parameters (`ParamsBiotSavart`) instead of just the
  long-range parameters (`ParamsLongRange`).

## [0.27.8] - 2025-06-12

### Changed

- Makie extension (for plotting): adapt to changes in Makie 0.23. We now require this version of Makie.

## [0.27.7] - 2025-06-10

### Changed

- `MinimalEnergy` mode: throw error if forcing is enabled. Also, store self-induced velocity
  (`vs`) in addition to actual filament velocity (`vL`).

### Added

- Add `Diagnostics.energy_injection_rate` for estimation of energy injection (or
  dissipation) rate due to external velocities / forcing.

- Add `Diagnostics.helicity_spectrum` for estimation of the helicity spectrum.

- Coarse-grained vorticity: print warning if smoothing scale is too small (compared to Ewald
  splitting scale) to guarantee good accuracy.

- Coarse-grained vorticity: allow computation from an existent _vorticity_ field (as opposed
  to a velocity field, which was already possible).

## [0.27.6] - 2025-06-03

### Added

- Add `BiotSavart.to_coarse_grained_vorticity!` and document how to interpolate results onto
  vortex positions.

### Changed

- `LongRangeCacheState` now stores the smoothing scale `σ` (the width of the Gaussian
  filter) instead of just a `smoothed` boolean stating whether the stored field has been
  smoothed.

## [0.27.5] - 2025-05-26

### Changed

- `FourierBandForcingBS`: properly estimate energy injection rate when forcing more than one
  mode.

- `Timestepping`: don't show default callback and affect functions in timings (`iter.to`) to
  reduce visual noise.

### Added

- Add `NoForcing` type. This may be used in `Timestepping.init` for explicitly disabling
  forcing (it's also the default).

## [0.27.4] - 2025-05-23

### Changed

- `FourierBandForcingBS`: the `ε_target` parameter now includes energy injected at
  wavevector $-\boldsymbol{k}$ due to Hermitian symmetry.

## [0.27.3] - 2025-05-23

### Fixed

- Fix `FourierBandForcing` with a time-varying normal fluid velocity and when running on GPUs.

### Added

- Add Biot–Savart based forcing method (`FourierBandForcingBS`), which allows to force
  specific wavevectors and depends only on the current vortex configuration. In particular,
  compared to other forcing types, it doesn't require imposing an external "normal" fluid
  velocity.

## [0.27.2] - 2025-05-14

### Added

- Add `affect_t!` option to `Timestepping.init`.
  This is similar to `affect!`, but allows to modify the state of the solver in a
  fine-grained manner, e.g. before every Runge–Kutta substep.
  As opposed to `affect!` and `callback`, an `affect_t!` function takes two arguments:
  an `iter::VortexFilamentSolver` and the current time `t::Real`.

- Add `mode = MinimalEnergy()` option to `Timestepping.init`.
  This evolves the vortices in "pseudo-time" such that the initial condition tends to a
  minimal energy configuration.
  For closed vortices, this simply makes the vortices shrink over time until they disappear.

## [0.27.1] - 2025-04-28

### Added

- Write current VortexPasta and Julia versions to checkpoint files.

## [0.27.0] - 2025-04-25

### Added

- Add `save_checkpoint` and `load_checkpoint` to `Timestepping` submodule.
  The former should be preferred to `FilamentIO.write_vtkhdf` when writing simulation results.

## [0.26.8] - 2025-04-08

### Added

- Refinement: add option to remove points in high curvature regions when using the `RefineBasedOnSegmentLength` criterion.

## [0.26.7] - 2025-03-07

### Added

- Store slip velocity (Fourier-filtered `vn - vs`) evaluated on filament nodes when using `FourierBandForcing`.
  This is now available as the `iter.v_ns` field of `VortexFilamentSolver` when this kind of forcing is used.

- More informative (and nicer looking) output of `println(iter::VortexFilamentSolver)`.

## [0.26.6] - 2025-03-05

### Fixed

- Fix possible issue when using `FourierBandForcing` on GPUs.

## [0.26.5] - 2025-03-05

### Added

- Implement `FourierBandForcing` as an alternative to `NormalFluidForcing`.
  This new forcing should be more localised in scale space, and therefore may
  be more relevant for injecting energy in turbulence simulations.

## [0.26.4] - 2025-03-04

### Changed

- In `BiotSavart.compute_on_nodes!`, the correction to the long-range
  self-interaction is now included in the `shortrange` computation instead of the `longrange` one.
  This might affect computations using schemes which split these both components (not sure these are very useful).

### Added

- Add `BiotSavart.get_longrange_field_fourier` function, which allows to obtain the
  latest computed field (vorticity, velocity, streamfunction) in Fourier space.

- Add basic text format to `FilamentIO`, as an alternative to the VTKHDF format.

## [0.26.3] - 2025-02-18

### Fixed

- Fix computation of energy spectrum on GPUs.
  The issue is probably linked to a recent update in KernelAbstractions.jl.

## [0.26.2] - 2025-02-17

### Added

- Allow NonuniformFFTs 0.7.

## [0.26.1] - 2025-01-27

### Changed

- Allow MakieCore 0.9.

## [0.26.0] - 2024-12-21

### Added

- Add option to inject filaments during a simulation. This should be done by
  calling the new `Timestepping.inject_filament!` function from within an `affect!` callback
  (see docs for details).

### Changed

- Reduce variability of results due to recent parallelisation of reconnections.
  Previously, the order of reconnections was somewhat arbitrary, depending on
  which threads tested which segment pairs and in which order. We now give
  priority to pairs which have the smallest distance.
  This means that results might slightly change compared to simulations with
  previous versions, since the priority given to reconnection candidates has
  changed (things are more deterministic now).

## [0.25.4] - 2024-12-19

### Changed

- Try to always parallelise (using threads) over filament nodes as opposed to separate filaments.
  This will speed up things in cases where there are few filaments, or even
  a single large filament along with many small ones.
  Also, we now parallelise reconnections (candidate finding) as well.

## [0.25.3] - 2024-12-16

### Fixed

- Workaround CUDA failure when transferring data between host and device.
  This could happen in simulations where the number of vortex points tended to increase over time
  (e.g. forced simulations).
  Seems to be related to an issue with memory-pinned CPU arrays whose size (and
  thus their location in memory) can change over time.
  See [CUDA.jl issue](https://github.com/JuliaGPU/CUDA.jl/issues/2594) for details.

## [0.25.2] - 2024-12-16

### Fixed

- Workaround possible precompilation failure when setting multiple CPU targets
  via the `JULIA_CPU_TARGET` environment variable.
  This is useful on HPC clusters where different nodes may have different CPU types.
  In particular, the problem was observed when setting `JULIA_CPU_TARGET=generic;sapphirerapids,clone_all;cascadelake,clone_all` on the Jean-Zay IDRIS cluster.

## [0.25.1] - 2024-12-09

### Changed

- **Breaking**: When a normal fluid is applied, the `vs` field of
  a `VortexFilamentSolver` (usually `iter.vs`) now represents the self-induced
  vortex velocity, while the `vL` field is the actual filament velocity
  including the contribution from mutual friction.

- Short-range parallelisation using threads: use finer-grained parallelisation
  at the level of the nodes of each filament. This is fine since the cost of
  computing all short-range interactions on a filament node is relatively large,
  so that the overhead of threading is relatively very small. Moreover, it's
  much faster when we have a small number of filaments (which, counterintuitively,
  can be relevant for turbulent cases). Also, we now avoid using `@threads :static`.

### Added

- Store self-induced velocity ($v_s$) in addition to the actual filament
  velocity ($v_L$) in `Timestepping.VortexFilamentSolver`.

- Allow saving state of `FourierBandVectorField` to HDF5 file and loading it back.
  This is useful in particular for restarts.

## [0.25.0] - 2024-12-06

### Added

- Add `VortexPasta.Forcing` module allowing to specify "forcing" methods.
  Currently, this includes the `NormalFluidForcing` type for modifying the
  vortex velocities via mutual friction with a (synthetic) normal fluid.
  Moreover, the `init` function (in `Timestepping`) now accepts a `forcing`
  argument which can be used in simulations.

- Add `VortexPasta.SyntheticFields` module allowing to create and evaluate
  synthetic vector fields (such as a normal fluid velocity field).

### Fixed

- Fix possible issue in reconnections when using velocity information (added in v0.24.8).
  When a reconnection was performed and a filament removed or added, the
  filament indices were not modified, which could lead to wrong velocities or
  even out-of-bounds accesses. Now, we find all reconnection pairs before doing
  any actual reconnections, which avoids this issue and simplifies a few things.

- Fix potential performance issues when using `Float32` on AVX512-capable CPUs.
  The issue affected short-range Biot–Savart computations, and more specifically
  the part where explicit SIMD vectorisation is used (for computing `erfc`
  among other things). We used to splat a `MMatrix` (from StaticArrays.jl)
  whose size could exceed 32 elements when using `Float32` on an AVX512-capable
  CPU, which seems to cause inference issues and bad performance.

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
