"""
    helicity_spectrum(iter::VortexFilamentSolver) -> (ks, Hk)
    helicity_spectrum(cache) -> (ks, Hk)

Compute helicity spectrum associated to vortex filament state.

Returns a tuple of vectors `(ks, Hk)` where `ks` contains the probed wavenumbers and `Hk`
the helicity associated to each wavenumber.

See also [`helicity_spectrum!`](@ref) for a non-allocating variant and for more details.
"""
function helicity_spectrum end

@doc raw"""
    helicity_spectrum!(Hk::AbstractVector, ks::AbstractVector, cache)

Compute helicity spectrum associated to vortex filament state.

For a single wave vector ``\bm{k}``, the associated helicity is defined as

```math
H(\bm{k}) ≡ \frac{V}{Δk} \bm{v}^*_{\bm{k}} ⋅ \bm{ω}_{\bm{k}}
= 2\frac{V}{Δk} \bm{k} ⋅ \left[ \bm{v}_{\text{r}} × \bm{v}_{\text{i}} \right].
```

In the last expression, the subscripts ``\text{r}`` and ``\text{i}`` respectively denote
real and imaginary parts.
That result can be easily found by doing some algebra, after replacing
``\bm{ω}_{\bm{k}} = \mathrm{i} \bm{k} × \bm{u}_{\bm{k}}``.
It explicitly shows that ``H(\bm{k})`` is purely real.

As usual, this function computes the sum of ``H(\bm{k})`` over spherical shells of radius
``k_n ≤ |\bm{k}| < k_n + Δk = k_{n + 1}``.

Defined as above, the total helicity in a (cubic) periodic domain of volume ``V = L^3`` is
``H = ∑_{\bm{k}} H(\bm{k}) Δk`` where ``Δk = 2π/L`` is the resolution in wavenumber space.

See [`energy_spectrum!`](@ref) for details on accepted `cache` arguments.
"""
function helicity_spectrum! end

function helicity_spectrum(cache; kws...)
    ks, Hk = init_spectrum(cache)
    helicity_spectrum!(Hk, ks, cache; kws...)
end

function helicity_spectrum!(
        Hk::AbstractVector, ks::AbstractVector, cache::LongRangeCache;
    )
    (; state,) = cache.common
    p = BiotSavart.get_parameters(cache)
    V = prod(p.Ls)  # domain volume
    from_velocity = state.quantity == :velocity && state.smoothing_scale == 0
    from_vorticity = state.quantity == :vorticity && state.smoothing_scale == 0
    if from_velocity
        _compute_spectrum!(Hk, ks, cache) do u⃗, k⃗, k², I
            # Note: û ⋅ ω̂ is expected to be purely real.
            # Also, the factor 2 is because _compute_spectrum applies a factor 0.5 (for
            # energy), which is not present in the helicity.
            ω⃗ = im * (k⃗ × u⃗)
            2V * real(u⃗ ⋅ ω⃗)
        end
    elseif from_vorticity
        _compute_spectrum!(Hk, ks, cache) do ω⃗, k⃗, k², I
            local β = ifelse(iszero(k²), one(k²), 1 / k²)
            local u⃗_lap = im * (k⃗ × ω⃗)  # this is k² * û
            2V * β * real(u⃗_lap ⋅ ω⃗)
        end
    else
        error(lazy"the current state of the long-range cache ($state) is currently not supported")
    end
    ks, Hk
end

helicity_spectrum!(Hk::AbstractVector, ks::AbstractVector, cache; kws...) =
    helicity_spectrum!(Hk, ks, get_long_range_cache(cache); kws...)
