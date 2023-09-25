# [The vortex filament model](@id VFM)

The standard **vortex filament model** (VFM) describes the motion of thin vortex lines in three-dimensional space.

Vortex lines are assumed to be very thin with respect to the scales of interest, such that they can be effectively described as spatial curves.

## [The Biot--Savart law](@id BiotSavart)

In the VFM, each vortex line induces a velocity field about it given by the [Biot--Savart law](https://en.wikipedia.org/wiki/Biot%E2%80%93Savart_law):

```math
\bm{v}(\bm{x}) =
\frac{Γ}{4π} ∮_{\mathcal{C}} \frac{(\bm{s} - \bm{x}) \times \mathrm{d}\bm{s}}{|\bm{s} - \bm{x}|^3}
```

where ``Γ`` is the vortex circulation (equal to ``κ`` for quantum vortices),
and ``\bm{s}`` is a point along the vortex line.
Here ``\mathcal{C}`` denotes the whole set of vortex lines in the system.

Mathematically, the above equation derives from the vorticity field:

```math
\bm{ω}(\bm{x}) ≡ \bm{\nabla} × \bm{v}(\bm{x})
= Γ ∮_{\mathcal{C}} δ(\bm{s} - \bm{x}) \, \mathrm{d}\bm{s}
```

where ``δ`` is Dirac delta function.
That is, the vorticity field is singular and localised at the locations of quantum vortices.

The Biot--Savart law describes in particular the motion induced by vortex filaments on themselves and on surrounding vortex lines.
The VFM thus describes the collective motion of a set of mutually-interacting vortex filaments which obey the Biot--Savart law.
Note that the Biot--Savart integral is singular when evaluated at a vortex location ``\bm{s}' ∈ \mathcal{C}``, and the integral must be desingularised by taking into account the finite thickness of the vortex core.
The VFM also accounts for *vortex reconnections*, which occur when two vortex segments are sufficiently close to each other and which change the topology of the vortex system.

## Desingularisation

In the VFM, one usually wants to evaluate the Biot--Savart law at locations ``\bm{x} = \bm{s}_0`` on the vortex.
It is clear that the Biot--Savart integral, as written above, is singular when evaluated at a point ``\bm{s}_0`` on the curve.

The divergence of the Biot--Savart integral is of course unphysical, and is related to the fact that the actual thickness of vortex lines not really infinitesimal but finite.
The standard way of accounting for the radius ``a`` of the vortex core is to split the integral into local and non-local parts:

```math
\bm{v}(\bm{s}_0) =
\frac{Γ}{4π} ∫_{\mathcal{C}_0} \frac{(\bm{s} - \bm{s}_0) \times \mathrm{d}\bm{s}}{|\bm{s} - \bm{s_0}|^3}
+
\frac{Γ}{4π} ∫_{\mathcal{C} ∖ \mathcal{C}_0} \frac{(\bm{s} - \bm{s}_0) \times \mathrm{d}\bm{s}}{|\bm{s} - \bm{s_0}|^3}
= \bm{v}_{\text{local}}(\bm{s}_0) + \bm{v}_{\text{non-local}}(\bm{s}_0)
```

Here ``\mathcal{C}_0`` denotes a portion of the set of curves ``\mathcal{C}`` which is in the neighbourhood of the point of interest ``\bm{s}_0``.

To illustrate this, the figure below shows a [trefoil knot](https://en.wikipedia.org/wiki/Trefoil_knot) curve, or rather its projection on the XY plane.
Note that here we represent a *discretised* version of the curve, where the number of degrees of freedom is finite and controlled by the positions of the markers (see [Generate figures](@ref) below for the code used to generate this figure).

![](trefoil_local.svg)

Here, to evaluate the velocity induced by the trefoil vortex on its point ``\bm{s}_i``, we split the curve into a local part ``\mathcal{C}_i`` (orange) and a non-local part ``\bar{\mathcal{C}}_i = \mathcal{C} ∖ \mathcal{C}_i`` (grey).
The non-local part is far from the singularity, so there is no need to modify the Biot--Savart integral as written above.
As for the local part, we can approximate it using a Taylor expansion of the Biot--Savart integral about ``\bm{s}_i`` and truncating the integral at a small distance ``ϵ ∝ a`` from the singularity.

More explicitly, from a Taylor expansion of ``\bm{s}`` close to ``\bm{s}_i = \bm{s}(ξ_i)``,
one can show that the Biot--Savart integrand has the Taylor expansion

```math
\frac{[\bm{s}(ξ) - \bm{s}_i] × \bm{s}'(ξ)}{|\bm{s}(ξ) - \bm{s}_i|^3}
≈ \frac{\bm{s}_i' × \bm{s}_i''}{2 (ξ - ξ_i)},
```

where ``ξ`` is the curve arc length and derivatives (primes) are with respect to ``ξ``.
Note that ``\bm{s}_i'`` and ``\bm{s}_i''`` are respectively the unit tangent and curvature vectors at ``\bm{s}_i`` (see figure).

Now, if one integrates e.g. from ``\bm{s}(ξ_i + ϵ)`` to ``\bm{s}_{i + 1} = \bm{s}(ξ_{i + 1})``, one gets:

```math
\bm{v}_i^+
≈ \frac{Γ}{4π} \frac{\bm{s}_i' × \bm{s}_i''}{2}
∫_{ξ_i + ϵ}^{ξ_{i + 1}} \frac{\mathrm{d}ξ}{ξ - ξ_i}
= \frac{Γ}{4π} \frac{\bm{s}_i' × \bm{s}_i''}{2}
\ln \left( \frac{ℓ^+}{ϵ} \right),
```

where ``ℓ^+ = ξ_{i + 1} - ξ_i`` is the length of the curve segment ``\bm{s}_i → \bm{s}_{i + 1}`` (see figure).
Doing something similar for the segment ``\bm{s}_{i - 1} → \bm{s}(ξ_i - ϵ)``, we obtain the local velocity:

```math
\bm{v}_{\text{local}}(\bm{s}_i)
= \frac{Γ}{4π} \bm{s}_i' × \bm{s}_i'' \ln \frac{\sqrt{ℓ^- ℓ^+}}{ϵ}
= \frac{Γ}{4π} \bm{s}_i' × \bm{s}_i'' \left[ \ln \frac{2\sqrt{ℓ^- ℓ^+}}{a} - Δ \right],
```

where we have taken the cut-off distance to be ``ϵ = \frac{e^Δ a}{2}`` [Saffman1993; §11.2](@cite).
Here ``Δ`` is a coefficient which depends on the form of the vorticity profile within the vortex core (see the [vortex ring tutorial](@ref Computing-the-vortex-ring-velocity)).

## Generate figures

### Trefoil knot figure

```@example
using CairoMakie
CairoMakie.activate!(type = "svg", pt_per_unit = 1.0)
Makie.set_theme!()

using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot

trefoil = define_curve(TrefoilKnot())
N = 32
i = (N ÷ 2) + 3
refinement = 8
colours = Makie.wong_colors()
f = Filaments.init(trefoil, ClosedFilament, N, CubicSplineMethod())

fig = Figure(Text = (fontsize = 24,))
ax = Axis(fig[1, 1]; aspect = DataAspect(), xlabel = "x", ylabel = "y")
hidexdecorations!(ax; label = false, ticklabels = false, ticks = false)
hideydecorations!(ax; label = false, ticklabels = false, ticks = false)

# Plot complete filament in grey
let color = (:grey, 0.6)
    plot!(ax, f; refinement, color, linewidth = 1.5)
    text!(
        ax, f(i ÷ 2 + 1, 0.6);
        text = L"\bar{\mathcal{C}}_{\!i}", align = (:right, :bottom), color,
        fontsize = 28,
    )
end

# Plot point of interest and local segment around it
let color = colours[2], arrowcol = colours[3]
    # Plot nodes (i - 1):(i + 1)
    scatter!(ax, getindex.(Ref(f), (i - 1):(i + 1)); color)
    text!(
        ax, f[i - 1] + Vec3(0.08, 0.0, 0.0);
        text = L"\mathbf{s}_{i - 1}", align = (:left, :center), color,
    )
    text!(
        ax, f[i] + Vec3(0.08, 0.0, 0.0);
        text = L"\mathbf{s}_i", align = (:left, :center), color,
    )
    text!(
        ax, f[i + 1] + Vec3(-0.08, 0.08, 0.0);
        text = L"\mathbf{s}_{i + 1}", align = (:right, :center), color,
    )

    # Plot local segments
    let ζs = range(0, 1; length = refinement + 1)
        lines!(ax, f.(i - 1, ζs); color, linewidth = 3)
        lines!(ax, f.(i, ζs); color, linewidth = 3)
        text!(
            f(i - 1, 0.5);
            text = L"ℓ^{-}", align = (:right, :bottom), color,
        )
        text!(
            f(i, 0.6) + Vec3(-0.08, 0, 0);
            text = L"ℓ^{+}", align = (:right, :center), color,
        )
        text!(
            ax, f(i - 1, 0.7) + Vec3(0.2, 0.0, 0.0);
            text = L"\mathcal{C}_{\!i}", align = (:left, :top), color,
            fontsize = 28,
        )
    end

    # Plot tangent and curvature vectors (assuming this is a 2D plot...)
    s⃗ = f[i]
    t̂ = f[i, UnitTangent()]
    ρ⃗ = f[i, CurvatureVector()]
    arrows!(ax, [s⃗[1]], [s⃗[2]], [t̂[1]], [t̂[2]]; color = arrowcol)
    arrows!(ax, [s⃗[1]], [s⃗[2]], [ρ⃗[1]], [ρ⃗[2]]; color = arrowcol)
    text!(
        ax, s⃗ + t̂ + Vec3(0.05, 0.0, 0.0);
        text = L"\mathbf{s}′", align = (:left, :center), color = arrowcol,
    )
    text!(
        ax, s⃗ + ρ⃗ + Vec3(0.0, 0.04, 0.0);
        text = L"\mathbf{s}″", align = (:center, :bottom), color = arrowcol,
    )
end

save("trefoil_local.svg", fig)
nothing
```

## References

A classical reference for the VFM is [Schwarz1985](@citet), while a more modern review is given by [Hanninen2014](@citet).
