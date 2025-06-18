@doc raw"""
    SmallScaleDissipationBS <: AbstractForcing
    SmallScaleDissipationBS(; kdiss, α, ε_target)

Dissipation based on Biot–Savart energetics at large wavenumbers (small scales).

This dissipation method works similarly to [`FourierBandForcingBS`](@ref). The main
difference is that it acts at a different range of scales (wavenumbers `k > kdiss`), and
that positive values of `α` or `ε_target` lead to dissipation and not injection.

As in [`FourierBandForcingBS`](@ref), one should pass _either_ `α` or `ε_target` but never
both. They should be positive for energy dissipation (negative values lead to energy
injection at small scales, which may lead to instabilities!):

- `α` (`\alpha`) is a non-dimensional coefficient which directly sets the amplitude of the dissipation term;

- `ε_target` (`\varepsilon_target`) has the units of an energy dissipation rate. In this case,
  the amplitude ``α`` will be adjusted over time in order to keep a roughly constant energy
  dissipation rate (which in general will _not_ be equal to `ε_target`, as discussed in
  [`FourierBandForcingBS`](@ref)).

# Extended help

## Dissipation term definition

The idea is to apply a "dissipation velocity" ``\bm{v}_{\text{diss}}`` to the vortices which
will extract energy from the small scales. Similarly to the forcing velocity in
[`FourierBandForcingBS`](@ref), such dissipation velocity can be written as:

```math
\bm{v}_{\text{diss}}(\bm{s}) = -α \bm{s}' × \bm{v}^{>}(\bm{s})
```

where ``\bm{v}^{>}(\bm{s})`` denotes the high-pass filtered superfluid velocity (from
Biot–Savart's law) at scales ``k > k_{\text{diss}}``.

In practice, this velocity is obtained as ``\bm{v}^{>} = \bm{v} - \bm{v}^{<}``, where
``\bm{v}`` and ``\bm{v}^{<}`` are the total and low-pass filtered superfluid velocities (up to
scale ``k_{\text{diss}}`` included).
This means that the Fourier-space fields must be sufficiently well resolved, up to the
chosen dissipative scale.
Evaluating ``\bm{v}^{<}`` at the vortex positions requires an extra type-2 NUFFT
(interpolation from Fourier to physical space), which is a bit costly.

## Energy dissipation rate

At scales ``k > k_{\text{diss}}``, the dissipation rate associated to this term is:

```math
ε_{\text{diss}} = -α \frac{Γ}{V} ∮ |\bm{s}' × \bm{v}^{>}|^2 \, \mathrm{d}ξ
```
"""
struct SmallScaleDissipationBS <: AbstractForcing
    # TODO
end
