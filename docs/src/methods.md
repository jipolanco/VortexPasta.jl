# Numerical methods

This page summarises some numerical aspects which should be taken into account when choosing solver parameters.

## Ewald summation for Biot--Savart

The main originality of the VortexPasta solver is that it adapts the [Ewald summation](https://en.wikipedia.org/wiki/Ewald_summation) method to accelerate the computation of the Biot--Savart law along vortex filaments.
See for example [Arnold2005](@citet) for a nice introduction to Ewald methods applied to the electrostatic interaction between point charges.

The basic idea of the method is to split the Biot--Savart integral into short- and long-range parts:
```math
\begin{align*}
    \bm{v}(\bm{x})
    &= \frac{Γ}{4π} ∮_{\mathcal{C}}
    \frac{(\bm{s} - \bm{x}) \times \mathrm{d}\bm{s}}{|\bm{s} - \bm{x}|^3}
    \\
    &= \frac{Γ}{4π} ∮_{\mathcal{C}}
    \left[ \textcolor{#D62728}{g^{<}(|\bm{s} - \bm{x}|)} + \textcolor{#1F77B4}{g^{>}(|\bm{s} - \bm{x}|)} \right]
    \frac{(\bm{s} - \bm{x}) \times \mathrm{d}\bm{s}}{|\bm{s} - \bm{x}|^3}
    \\
    &= \textcolor{#D62728}{\bm{v}^{<}(\bm{x})} +
       \textcolor{#1F77B4}{\bm{v}^{>}(\bm{x})}
\end{align*}
```
where the *short-range* and *long-range* scalar functions ``g^<(r)`` and ``g^>(r)`` have the properties:
1. ``g^<(r) + g^>(r) = 1`` for all ``r > 0``;
2. ``g^<(r)`` decays exponentially with ``r``, so that long-range interactions can be neglected beyond some cut-off distance ``r_\text{cut}`` when computing the short-range velocity ``\bm{v}^<``;
3. ``g^>(r) / r^2`` is non-singular and smooth at ``r = 0``, which implies that ``g^>(r)`` must *quickly* tend to 0 as ``r → 0``. In periodic domains, this enables the use of the fast Fourier transform (FFT) to efficiently estimate long-range interactions.

One can show that the following pair of functions obeys all these properties:
```math
\begin{align*}
    g^<(r) &= \operatorname{erfc}(αr) + \frac{2αr}{\sqrt{π}} e^{-α^2 r^2}
    \\
    g^>(r) &= \operatorname{erf}(αr) - \frac{2αr}{\sqrt{π}} e^{-α^2 r^2}
\end{align*}
```
Here ``\operatorname{erf}`` and ``\operatorname{erfc}`` are respectively the [error function](https://en.wikipedia.org/wiki/Error_function) and the [complementary error function](https://en.wikipedia.org/wiki/Error_function#Complementary_error_function), which satisfy ``\operatorname{erf}(z) + \operatorname{erfc}(z) = 1``.
Above we have introduced the **Ewald splitting parameter** ``α``, which is the inverse of a length scale.
This parameter is purely numerical and, in theory, its choice has no impact on the final result ``\bm{v}(\bm{x})``.
In practice this is the case as long as other numerical parameters are well chosen (see discussion further below).

The two splitting functions are plotted below.
In the horizontal axis, the scale ``r`` is non-dimensionalised by the splitting parameter ``α``.

```@example
using SpecialFunctions: erf, erfc
using CairoMakie
CairoMakie.activate!(type = "svg", pt_per_unit = 1.0)  # hide
Makie.set_theme!()  # hide
rs = 2.0.^(-4:0.1:3)
gs(αr) = erfc(αr) + 2αr / sqrt(π) * exp(-αr^2)  # short-range
gl(αr) =  erf(αr) - 2αr / sqrt(π) * exp(-αr^2)  # long-range
xticks = LogTicks(-4:1:3)
yticks = LogTicks(-16:4:0)
fig = Figure(resolution = (600, 400), fontsize = 18)
ax = Axis(fig[1, 1]; xticks, yticks, xscale = log2, yscale = log10, xlabel = L"αr", ylabel = "Splitting function")
ylims!(ax, 1e-17, 4)
lines!(ax, rs, gs.(rs); label = L"g^<(r)")
lines!(ax, rs, gl.(rs); label = L"g^>(r)")
let rs = 2.0.^(range(-4, -1; length = 3)), color = :grey20  # plot ~r^3 slope
    ys = @. 0.2 * rs^3
    lines!(ax, rs, ys; linestyle = :dash, color)
    text!(ax, rs[2], ys[2]; text = L"r^3", align = (:left, :top), color)
end
axislegend(ax; position = (0, 0), labelsize = 20)
fig
```

For small ``αr``, the long-range splitting function goes to zero as ``r^3``, consistently with its Taylor expansion, ``g^>(r) = \frac{4}{3 \sqrt{π}} (αr)^3 + \mathcal{O}(r^5)``.
This means that, as we wanted, ``g^>(r) / r^2`` is non-singular and smooth at ``r = 0``.

As for the short-range splitting function, it is dominant for small ``αr``, while it decays exponentially to 0 for large ``αr``.
In particular, for ``r ≳ 6/α``, its value decays below about ``10^{-15}``, meaning that it is safe to set the cut-off distance ``r_\text{cut}`` around this value.

### Computation of long-range velocity

The long-range velocity ``\bm{v}^>`` has a simple physical interpretation.
Indeed, it can be shown, by differentiating the splitting functions defined above,
that its associated long-range vorticity ``\bm{ω}^>(\bm{x}) ≡ \bm{∇} × \bm{v}^>(\bm{x})``
is nothing else that a Gaussian-filtered version of the actual vorticity field induced by the vortex filaments (which is singular).

More precisely, the long-range vorticity is given by the convolution
``\bm{ω}^> = H ∗ \bm{ω}``, where ``H(\bm{r}) = (α / \sqrt{π})^3 \, e^{-α^2 r^2}``
is a 3D Gaussian kernel.
Intuitively, this means that the long-range velocity ``\bm{v}^>`` corresponds to the velocity induced by a coarse-grained version of the vortex filaments.
In this view, vortices are not "infinitesimal" anymore, but are closer to what we are used to see in classical viscous fluids.
In periodic domains, such a smooth vorticity field can be accurately and efficiently expanded in Fourier series, and the curl operator can be readily inverted in Fourier space to obtain the coarse-grained velocity field ``\bm{v}^>``.

#### 1. Estimating the vorticity in Fourier space

To do this, we first expand the (*actual*) vorticity field in Fourier series:

```math
\bm{ω}(\bm{x}) = Γ ∮_{\mathcal{C}} δ(\bm{x} - \bm{s}) \, \mathrm{d}\bm{s}
= ∑_{\bm{k}} \hat{\bm{ω}}(\bm{k}) \, e^{i \bm{k} ⋅ \bm{x}}
\quad\text{for } \bm{k} = \frac{2π \bm{n}}{L}, \bm{n} ∈ \mathbb{Z}^3
```

Here ``L`` is the domain period (assumed the same in all directions for simplicity).
The Fourier coefficients ``\hat{\bm{ω}}(\bm{k})`` are given by

```math
\hat{\bm{ω}}(\bm{k})
= \frac{1}{L^3} ∫ \bm{ω}(\bm{x}) \, e^{-i \bm{k} ⋅ \bm{x}} \, \mathrm{d}\bm{x}
= \frac{Γ}{L^3} ∮_{\mathcal{C}} e^{-i \bm{k} ⋅ \bm{s}} \, \mathrm{d}\bm{s}
```

To estimate this integral, we discretise the curves defining the support ``\mathcal{C}`` of the vortex filaments, leading to a weighted sum over discrete points ``\bm{x}_j``:

```math
\hat{\bm{ω}}(\bm{k})
≈ \frac{Γ}{L^3} ∑_{j} \bm{w}_j \, e^{-i \bm{k} ⋅ \bm{x}_j}
```

The weights are related to the length of the discrete segments and to the quadrature rule used to estimate the integrals over segments (more details soon...).
The important point here is that, since the discrete points are not generally on an equispaced grid, one cannot directly use the fast Fourier transform (FFT) to efficiently obtain these coefficients.
Nevertheless, they can be efficiently and accurately estimated using the [non-uniform FFT](https://en.wikipedia.org/wiki/Non-uniform_discrete_Fourier_transform#Nonuniform_fast_Fourier_transform) (NUFFT) algorithm.
More precisely, this corresponds to a [type-1 NUFFT](https://finufft.readthedocs.io/en/latest/math.html), from non-uniform points in physical space to uniform wavenumbers ``\bm{k}`` in Fourier space.

#### 2. Coarse-grained vorticity and velocity in Fourier space

Once we have obtained the ``\hat{\bm{ω}}(\bm{k})`` coefficients, obtaining the Fourier coefficients for the coarse-grained vorticity and coarse-grained velocity is straightforward.
The former is obtained by convolution with a Gaussian, which is simply a product in Fourier space:

```math
\hat{\bm{ω}}^>(\bm{k})
= \hat{\bm{ω}}(\bm{k}) \, e^{-k^2 / 4α^2}, \quad\text{where } k = |\bm{k}|
```

Similarly, the curl operator can be easily inverted in Fourier space to get the coarse-grained velocity:

```math
\hat{\bm{v}}^>(\bm{k})
= \frac{i \bm{k}}{k^2} × \hat{\bm{ω}}^>(\bm{k})
= \frac{i \bm{k}}{k^2} × \hat{\bm{ω}}(\bm{k}) \, e^{-k^2 / 4α^2}
\quad\text{for } k ≠ 0
```

!!! warning "Zero mean vorticity condition"

    The velocity is well-defined only if ``\hat{\bm{ω}}(\bm{0}) = \bm{0}``, that is, if the mean vorticity of the vortex filament system is zero.
    Otherwise we get division by zero, which is related to the fact that we're dealing with an infinitely periodic system and energy diverges with a non-zero mean vorticity.

    This condition is automatically satisfied when dealing with closed vortex filaments.
    This may however not be the case for infinite filaments (for instance, putting a single straight infinite filament in the domain is ill-defined).

#### 3. Notes on required resolution

Above we have assumed that we can evaluate Fourier coefficients for any wavenumber ``\bm{k}``.
In fact, for practical reasons, we cannot evaluate all coefficients ``\hat{\bm{ω}}(\bm{k})`` for every possible ``\bm{k}``, and we need to set the number of wavenumbers ``N`` to compute in each direction (this is a parameter of NUFFT algorithms).
In other words, we need to truncate the estimations at some maximum wavenumber, namely the Nyquist frequency, which is related to ``N`` by ``k_\text{max} = π N / L``.

Similarly to the cut-off distance in physical space, one can expect that the appropriate value of ``k_\text{max}`` (which is an inverse length scale) to get good accuracy should be proportional to the Ewald splitting parameter ``α``.
A rule of thumb is to choose a wavenumber at which the Gaussian factor ``\exp \! \left(-k_{\text{max}}^2 / 4α^2 \right)`` matches the desired accuracy.
For instance, at ``k_{\text{max}} = 8α``, this factor has dropped to about ``10^{-7}``.
Of course, one can vary the ``k_\text{max} / α`` ratio depending on the wanted accuracy.

#### 4. Physical velocity at filament locations

The last step is to evaluate, from the coarse-grained velocity ``\hat{\bm{v}}^>(\bm{k})`` in Fourier space, the *physical* coarse-grained velocity ``\bm{v}^>(\bm{s})`` on vortex filament locations ``\bm{s}`` (which, once again, are generally not on a regular grid).
This operation can be written as:

```math
\bm{v}^>(\bm{s}_j) = ∑_{\bm{k}} \hat{\bm{v}}^>(\bm{k}) \, e^{i \bm{k} ⋅ \bm{s}_j}
```

for a set of locations ``\bm{s}_j ∈ \mathcal{C}``.
Note that in practice this sum is truncated to the chosen ``k_{\text{max}}``.

This operation can be efficiently computed using a type-2 NUFFT (from uniform wavenumbers ``\bm{k}`` to non-uniform locations ``\bm{x}``), which can be understood as an interpolation on the chosen points.

## Estimating line integrals
