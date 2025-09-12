export energy_flux

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
function flux_wavenumbers(cache::LongRangeCache, Nk::Integer)
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
    ks = flux_wavenumbers(cache, Nk)
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

        for j in eachindex(vs_filtered, fs)
            Filaments.update_coefficients!(vs_filtered[j]; knots = Filaments.knots(fs[j]))
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
