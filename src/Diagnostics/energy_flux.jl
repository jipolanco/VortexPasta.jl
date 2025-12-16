export energy_flux, energy_transfer_matrix

using Adapt: Adapt, adapt

@doc raw"""
    energy_flux(iter::VortexFilamentSolver, Nk::Integer; quad = nothing) -> (ks, fluxes)
    energy_flux(iter::VortexFilamentSolver, ks::AbstractVector; quad = nothing) -> (ks, fluxes)

Compute the energy "flux" across scales due to various velocity terms.

The second argument can be either the wanted number of output wavenumbers `Nk`, or the
actual output wavenumbers `ks`. In the first case, this function will choose
logarithmically-spaced wavenumbers `ks` across which the fluxes will be computed.
In general, the length of `ks` will be slightly lower than `Nk` to avoid possible
repetitions.

Note that the cost of computing the fluxes across each wavenumber is very high
(we basically perform one NUFFT per wavenumber), so `Nk` should not be too large,
and this function should not be called too often.

This function returns:

- `ks::AbstractVector`: a vector of wavenumbers;

- `fluxes::NamedTuple`: a list of energy fluxes across scales due to different velocity terms.
  Each element is a vector of real values.
  For example, `fluxes.vs` is a vector containing the energy flux across scales due to the
  Biot–Savart velocity (see Extended help below).

# Extended help

Depending on the fields available in `iter`, the returned `fluxes` can include the fields:

- `vs`: energy flux across scale ``k`` due to the Biot–Savart self-induced velocity;

- `vinf`: estimated non-local energy flux from scales ``≤ k`` to scales ``k ∼ ∞``;

- `vf`: accumulated energy injection rate (by external forcing term) up to scale ``k``;

- `vdiss`: accumulated energy dissipation rate (by external dissipation term) up to scale ``k``.

Note that all of these fluxes are properly signed according to the usual convention (i.e. they are ready to plot).
In particular, both `vf` and `vdiss` are positive if the respective terms inject and
dissipate energy as expected, while `vs` is positive in the presence of a direct energy cascade.

## Definitions

The energy flux (`fluxes.vs`) across scale ``k`` is defined as:

```math
Π_k = - \frac{Γ}{V} ∮ \left[ \bm{s}' × \bm{v}_{\text{s}}^{<} \right] ⋅ \bm{v}_{\text{s}} \, \mathrm{d}ξ,
```

where ``\bm{v}_{\text{s}}`` is the self-induced vortex velocity (according to the Biot–Savart law)
and the superscript ``<`` denotes a low-pass filter up to scale ``k``.
Note that this looks very similar to the [`energy_injection_rate`](@ref) definition, and
this is because ``Π_k`` can indeed be interpreted as the dissipation of large-scale energy
(associated to ``\bm{v}_{\text{s}}^{<}``) by the full Biot–Savart velocity ``\bm{v}_{\text{s}}``.

Similarly, the estimated non-local energy flux (`fluxes.vinf`) from scales ``≤ k`` to scales ``k ∼ ∞``
is defined as:

```math
\begin{align}
Π_k^{∞}
&= - \frac{Γ}{V} ∮ \left[ \bm{s}' × \bm{v}_{\text{s}}^{<} \right] ⋅ \bm{v}_{\text{s}}^{∞} \, \mathrm{d}ξ
\\
&= - \frac{Γ^2}{4π V} ∮ \bm{s}'' ⋅ \bm{v}_{\text{s}}^{<} \, \mathrm{d}ξ,
\end{align}
```

where ``\bm{v}_{\text{s}}^{∞} ≡ \frac{Γ}{4π} \bm{s}' × \bm{s}''`` is an estimate for the
"singular" part of the self-induced velocity, associated to very high wavenumbers ``k``.
This is only an estimate which allows to have an idea of the non-local energy transfers from
large to small scales. Also note the close similarity of this expression with the filament
stretching rate (see [`stretching_rate`](@ref)). Indeed, ``Π_k^{∞}`` is proportional to the
increase of filament length due to the large-scale velocity ``\bm{v}_{\text{s}}^{<}``.

If an external dissipation term is included, the accumulated energy dissipation rate (`fluxes.vdiss`) is given by:

```math
D_k = - \frac{Γ}{V} ∮ \left[ \bm{s}' × \bm{v}_{\text{s}}^{<} \right] ⋅ \bm{v}_{\text{diss}} \, \mathrm{d}ξ,
```

where ``\bm{v}_{\text{diss}}`` is the dissipative velocity.

Finally, the accumulated energy injection rate (`fluxes.vf`) is (note the ``+`` sign!):

```math
F_k = + \frac{Γ}{V} ∮ \left[ \bm{s}' × \bm{v}_{\text{s}}^{<} \right] ⋅ \bm{v}_{\text{f}} \, \mathrm{d}ξ,
```

where ``\bm{v}_{\text{f}}`` is the forcing velocity.
"""
function energy_flux end

# Callable struct to be used as a callback.
struct DropLargeWavenumbers{T <: Real, WavenumberVectors <: NTuple} <: Function
    kcut  :: T
    kcut² :: T
    ks :: WavenumberVectors
end

# For some reason, this seems to be needed for compilation on CUDA.
# (Even though I think it doesn't really do anything and doesn't change any types, since
# nothing is really allocated in CPU or GPU memory.)
function Adapt.adapt_structure(to, f::DropLargeWavenumbers)
    DropLargeWavenumbers(
        f.kcut,
        f.kcut²,
        map(k -> adapt(to, k), f.ks),
    )
end

DropLargeWavenumbers(kcut, ks) = DropLargeWavenumbers(kcut, kcut^2, ks)

@inline function (f::DropLargeWavenumbers)(û::NTuple{N}, idx::NTuple{N}) where {N}
    (; ks, kcut²,) = f
    k⃗ = ntuple(d -> @inbounds(ks[d][idx[d]]), Val(N))
    k² = sum(abs2, k⃗)
    û_zero = ntuple(d -> @inbounds(zero(û[d])), Val(N))::typeof(û)
    if k² > kcut²
        û_zero
    else
        û
    end
end

# Determine wavenumbers over which the flux will be computed.
function flux_wavenumbers_logspaced(cache::LongRangeCache, Nk::Integer)
    ks_full = init_spectrum_wavenumbers(cache)  # = 0:Δk:kmax
    Δk = ks_full[begin + 1] - ks_full[begin]
    ks_log = collect(logrange(ks_full[begin + 1], ks_full[end]; length = Nk))  # logarithmically-spaced wavenumbers
    unique!(floor.(ks_log ./ Δk) .* Δk)  # remove repeated and "fractional" wavenumbers
end

function energy_flux(cache_in, fs, velocities::NamedTuple, Nk_or_ks, params::ParamsBiotSavart; kws...)
    cache = get_long_range_cache(cache_in)::LongRangeCache
    energy_flux(cache, fs, velocities, Nk_or_ks, params; kws...)
end

function energy_flux(cache::LongRangeCache, fs, velocities::NamedTuple, Nk::Integer, params::ParamsBiotSavart; kws...)
    ks = flux_wavenumbers_logspaced(cache, Nk)
    energy_flux(cache, fs, velocities, ks, params; kws...)
end

function energy_flux(
        cache::LongRangeCache, fs, velocities::NamedTuple,
        ks::AbstractVector, params::ParamsBiotSavart;
        quad = nothing, vs_buf = similar(first(velocities).field),
    )
    (; wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache)
    @assert state.quantity == :velocity
    @assert state.smoothing_scale == 0

    # Initialise flux vectors.
    fluxes = map(_ -> similar(ks), velocities)::NamedTuple

    # Iterate over target wavenumbers.
    vs_filtered = vs_buf
    for (i, k_flux) in pairs(ks)
        callback = DropLargeWavenumbers(k_flux, wavenumbers)
        BiotSavart.interpolate_to_physical!(callback, cache)
        BiotSavart.copy_long_range_output!(vs_filtered, cache)

        @sync for chunks in FilamentChunkIterator(fs; full_vectors = true)
            Threads.@spawn for (j, _, _) in chunks
                Filaments.update_coefficients!(vs_filtered[j]; knots = Filaments.knots(fs[j]))
            end
        end

        # Compute dissipation of large-scale energy by each input velocity
        for n in eachindex(values(velocities), values(fluxes))
            # The integrand is [s′ × vs] ⋅ vL == [s′ × vs^>] ⋅ velocities[n]
            vs = vs_filtered
            # In principle this is type unstable, since vL can either be a vector of velocities or CurvatureVector().
            # In practice this doesn't seem to be an issue as the loop is probably unrolled.
            vL = velocities[n].field
            s = velocities[n].sign  # -1 if dissipation, +1 if injection
            fluxes[n][i] = s * Diagnostics.energy_injection_rate(fs, vL, vs, params; quad)
        end
    end

    ks, fluxes
end

## ========================================================================================== ##

@doc raw"""
    energy_transfer_matrix(iter::VortexFilamentSolver, Nk::Integer; quad = nothing) -> (ks, transfers)
    energy_transfer_matrix(iter::VortexFilamentSolver, ks::AbstractVector; quad = nothing) -> (ks, transfers)

Compute energy transfers between different Fourier shells.

Returns a square antisymmetric matrix of dimensions `(Nk + 1, Nk + 1)` where `Nk` is the length of
the returned `ks` vector (which might be smaller than the input `Nk`, if there are not
enough resolved wavenumbers). Note that the `ks` vector does not include the mode `k = 0`.

This function defines `Nk + 1` shells, such that the `i`-th shell contains wavenumbers `ks[i - 1] < k ≤ ks[i]`,
with the conventions `ks[0] = 0` and `ks[Nk + 1] = ∞`. In other words, the last "shell"
contains all small wavenumbers `k > ks[end]`.

See [`energy_flux`](@ref) for more details on the accepted arguments.

## Definitions

The energy transfer between Fourier shells A and B is defined as

```math
T_{AB} = \frac{Γ}{V} ∮ \left[ \bm{s}' × \bm{v}_{\text{s}}^{\text{A}} \right] ⋅ \bm{v}_{\text{s}}^{\text{B}} \, \mathrm{d}ξ,
```

where ``\bm{v}_{\text{s}}^{\text{A}}`` and ``\bm{v}_{\text{s}}^{\text{B}}`` are Fourier
band-pass-filtered fields, including only the coefficients within the chosen shells.
"""
function energy_transfer_matrix end

# Callable struct to be used as a callback (similar to DropLargeWavenumbers).
struct KeepFourierShell{T <: Real, WavenumberVectors <: NTuple} <: Function
    kmin  :: T
    kmin² :: T
    kmax  :: T
    kmax² :: T
    ks :: WavenumberVectors
end

function Adapt.adapt_structure(to, f::KeepFourierShell)
    KeepFourierShell(
        f.kmin,
        f.kmin²,
        f.kmax,
        f.kmax²,
        map(k -> adapt(to, k), f.ks),
    )
end

KeepFourierShell(kmin, kmax, ks) = KeepFourierShell(kmin, kmin^2, kmax, kmax^2, ks)

@inline function (f::KeepFourierShell)(û::NTuple{N}, idx::NTuple{N}) where {N}
    (; ks, kmin², kmax²) = f
    k⃗ = ntuple(d -> @inbounds(ks[d][idx[d]]), Val(N))
    k² = sum(abs2, k⃗)
    û_zero = ntuple(d -> @inbounds(zero(û[d])), Val(N))::typeof(û)
    if kmin² < k² <= kmax²
        û
    else
        û_zero
    end
end

function energy_transfer_matrix(cache_in, fs, vs, args...; kws...)
    @assert !(cache_in isa LongRangeCache)
    cache = get_long_range_cache(cache_in)
    energy_transfer_matrix(cache, fs, vs, args...; kws...)
end

function energy_transfer_matrix(
        cache::LongRangeCache, fs::AbstractVector, vs::AbstractVector, Nk::Integer, params::ParamsBiotSavart;
        kws...
    )
    ks = flux_wavenumbers_logspaced(cache, Nk)  # logarithmically-spaced wavenumbers
    energy_transfer_matrix(cache, fs, vs, ks, params; kws...)
end

function energy_transfer_matrix(
        cache::LongRangeCache, fs::AbstractVector, vs::AbstractVector, ks::AbstractVector, params::ParamsBiotSavart;
        quad = nothing,
    )
    (; wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache)
    @assert state.quantity == :velocity
    @assert state.smoothing_scale == 0

    Nk = length(ks)
    Base.require_one_based_indexing(ks)

    # Initialise filament velocities for each wavenumber shell
    nshells = Nk + 1
    vs_shells = [similar(vs) for _ in 1:nshells]
    vs_large = last(vs_shells)  # accumulated large-scale velocity

    # Compute shell-filtered velocities
    k_prev = zero(eltype(ks))
    for (i, k_shell) in enumerate(ks)
        # Shell limits
        kmin = k_prev
        kmax = k_shell
        callback = KeepFourierShell(kmin, kmax, wavenumbers)
        BiotSavart.interpolate_to_physical!(callback, cache)
        vs_shell = vs_shells[i]
        BiotSavart.copy_long_range_output!(vs_shell, cache)
        if i == 1
            vs_large .= vs_shell
        else
            vs_large .+= vs_shell  # vs_large now includes the velocity for 0 < k < k_shell
        end
        k_prev = k_shell
    end

    vs_shells[end] .= vs .- vs_large  # the last shell now has the small-scale velocity (wavenumbers k > ks[end])

    # TODO: parallelise over shells or over filament chunks?
    Threads.@threads for vs_shell in vs_shells
        for n in eachindex(vs_shell, fs)
            Filaments.update_coefficients!(vs_shell[n]; knots = Filaments.knots(fs[n]))
        end
    end

    # Compute energy transfer matrix
    transfers = similar(ks, (nshells, nshells))
    Threads.@threads for j in eachindex(vs_shells)
        transfers[j, j] = 0  # by definition diagonal values are zero
        for i in (j + 1):lastindex(vs_shells)
            # Don't parallelise inside of energy_injection_rate, since we're already running in parallel.
            Tij = Diagnostics.energy_injection_rate(fs, vs_shells[j], vs_shells[i], params; quad, nthreads = 1)
            transfers[i, j] = Tij
            transfers[j, i] = -Tij  # matrix is antisymmetric
        end
    end

    ks, transfers
end
