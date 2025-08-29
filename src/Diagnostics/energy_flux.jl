export energy_flux

using Adapt: Adapt, adapt

@doc raw"""
    energy_flux(iter::VortexFilamentSolver, Nk::Integer; quad = nothing) -> (ks, fluxes)

Compute the energy "flux" across scales due to various velocity terms.

The `Nk` argument is the maximum number of wavenumbers that will be considered.
This function will then define a vector of logarithmically-spaced wavenumbers `ks`
across which the fluxes will be computed.
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

Similarly, the accumulated energy dissipation rate (`fluxes.vdiss`) is given by:

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

function energy_flux(cache_in, velocities::NamedTuple, Nk::Integer, params::ParamsBiotSavart; kws...)
    cache = get_long_range_cache(cache_in)::LongRangeCache
    energy_flux(cache, velocities, Nk, params; kws...)
end

function energy_flux(
        cache::LongRangeCache, velocities::NamedTuple, Nk::Integer, params::ParamsBiotSavart;
        quad = nothing, vs_buf = similar(first(velocities)),
    )
    (; wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache)
    @assert state.quantity == :velocity
    @assert state.smoothing_scale == 0
    
    # Determine target wavenumbers.
    ks_full = init_spectrum_wavenumbers(cache)  # = 0:Δk:kmax
    Δk = ks_full[begin + 1] - ks_full[begin]
    ks_log = collect(logrange(ks[begin + 1], ks[end]; length = Nk))  # logarithmically-spaced wavenumbers
    ks_flux = unique!(floor.(ks_log ./ Δk) .* Δk)  # remove repeated and "fractional" wavenumbers

    # Initialise flux vectors.
    fluxes = map(_ -> similar(ks_flux), velocities)

    # Iterate over target wavenumbers.
    vs_filtered = vs_buf
    for (i, k_flux) in pairs(ks_flux)
        callback = DropLargeWavenumbers(k_flux, wavenumbers)
        BiotSavart.interpolate_to_physical!(callback, cache)
        BiotSavart.copy_long_range_output!(vs_filtered, cache)
        for i in eachindex(vs_filtered, fs)
            Filaments.update_coefficients!(vs_filtered[i]; knots = Filaments.knots(fs[i]))
        end
        # Compute dissipation of large-scale energy by each input velocity
        for n in eachindex(velocities, fluxes)
            # The integrand is [s′ × vs] ⋅ vL == - [s′ × vs^>] ⋅ velocities[n]
            vL = vs_filtered
            vs = velocities[n]
            fluxes[n][i] = Diagnostics.energy_injection_rate(fs, vL, vs, params; quad)
        end
    end

    # Energy injection term has a negative sign by convention (since it injects energy).
    if hasproperty(fluxes, :vf)
        fluxes.vf .*= -1
    end

    ks, fluxes
end
