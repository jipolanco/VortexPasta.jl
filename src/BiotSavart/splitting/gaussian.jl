@doc raw"""
    GaussianSplitting{T <: AbstractFloat} <: AbstractEwaldSplitting
    GaussianSplitting(; Ls, β, α, rcut, Ns)

Gaussian splitting kernel for Ewald summation.

This is the standard splitting kernel traditionally used in Ewald methods.

# Parameters

- `Ls::NTuple{3, T}`: domain period in each Cartesian direction (**mandatory**),
- `β::T`: nondimensional accuracy coefficient,
- `α::T`: splitting parameter (an inverse length scale),
- `rcut::T`: cut-off distance for short-range computations,
- `Ns::NTuple{3, Int}`: size of FFT grid for long-range computations.

Note that one does _not_ need to pass all of these parameters (see below).

# Construction

There are several ways of constructing a `GaussianSplitting`.
Note that the period `Ls` is always required.

The **recommended way** is to pass the nondimensional accuracy coefficient `β` and _one of_
`α`, `rcut` or `Ns`. For example:

```jldoctest
julia> Ls = (2π, 2π, 2π);
julia> β = 3.5;  # for 6-digit accuracy
julia> Ns = (128, 128, 128);
julia> splitting = GaussianSplitting(; Ls, β, Ns)
GaussianSplitting{Float64, 3} with:
 - Domain period:               Ls = (6.283185307179586, 6.283185307179586, 6.283185307179586)
 - Ewald splitting parameter:   α  = 9.0
 - Short-range cut-off:         r_cut = 0.3888888888888889 (rcut/Lmin = 0.061893588980181526)
 - Long-range resolution:       Ns = (128, 128, 128) (kmax = 63.0)
 - Short-range accuracy coeff.: β_shortrange = 3.5
 - Long-range cut-off coeff.:   β_longrange = 3.5
```

As detailed in [Polanco2025](@citet), the splitting and truncation parameters will then be set such that
``r_{\text{cut}} = β / α`` and ``k_{\text{max}} = 2 α β``.
One can use the following table to choose `β` according to the wanted precision:

| Precision digits |  Ewald ``β`` |
| :--------------: | :----------: |
|         3        |     2.0      |
|         4        |     2.5      |
|         6        |     3.5      |
|         8        |     4.0      |
|        10        |     4.5      |
|        12        |     5.0      |
|        14        |     5.5      |

For more control, one can pass all of `α`, `rcut` and `Ns`, in which case `β` will be ignored.
But doing this is not recommended as it can lead to loss of accuracy.

# Extended help

## Kernel definitions

This kernel splits the 3D Green's function ``G(\bm{r}) = 1 / (4πr)`` into the near- and far-field
contributions:

```math
G(\bm{r}) = G^{\text{(n)}}(\bm{r}) + G^{\text{(f)}}(\bm{r}) =
\frac{\operatorname{erfc}(αr)}{4πr} + \frac{\operatorname{erf}(αr)}{4πr}
```

where ``\operatorname{erf}(x) = \frac{2}{\sqrt{\pi}} \int_0^x e^{-t^2} \, \mathrm{d}t`` is the
[error function](https://en.wikipedia.org/wiki/Error_function) and ``\operatorname{erfc}(x) = 1 - \erf(x)``.
Here ``α`` is Ewald's splitting parameter (an inverse length scale).

As a result, the Biot--Savart kernel $\bm{\nabla}G(\bm{r}) = -\bm{r} / (4πr^3)$ is split as
$\bm{\nabla}G(\bm{r}) = \bm{\nabla}G^{\text{(n)}}(\bm{r}) + \bm{\nabla}G^{\text{(f)}}(\bm{r})$ with:

```math
\begin{align*}
    \bm{\nabla}G^{\text{(n)}}(\bm{r})
    &= -\frac{\bm{r}}{4\pi r^3} \left[ \operatorname{erfc}(αr) + \frac{2αr}{\sqrt{π}} \, e^{-(α r)^2} \right],
    \\
    \bm{\nabla}G^{\text{(f)}}(\bm{r})
    &= -\frac{\bm{r}}{4\pi r^3} \left[ \operatorname{erf}(αr) - \frac{2αr}{\sqrt{π}} \, e^{-(α r)^2} \right].
\end{align*}
```

Far-ranged fields can be interpreted as those induced by a Gaussian-filtered vorticity field
since ``-\nabla^2 G^{\text{f}}(\bm{r}) = (\alpha / \sqrt{\pi})^3 \, e^{-(\alpha r)^2} ≡ φ(\bm{r})``.
The Fourier transform of this convolution kernel is given by

```math
\hat{\varphi}(\bm{k}) = e^{-k^2 / (4 \alpha^2)}
```

## Kernel truncation

Given the accuracy and splitting parameters ``β`` and ``α``, the required physical and Fourier space
truncations are ``r_{\text{cut}} = β / α`` and ``k_{\text{max}} = 2 α β`` to achieve the desired accuracy.
"""
struct GaussianSplitting{T <: AbstractFloat, N} <: AbstractEwaldSplitting
    Ls::NTuple{N, T}
    α::T
    rcut::T
    Ns::Dims{N}
end

function GaussianSplitting(; Ls::NTuple{N, T}, β = nothing, α = nothing, rcut = nothing, Ns = nothing) where {N, T}
    GaussianSplitting(Ls, _gaussian_splitting_params(Ls, β, α, rcut, Ns)...)
end

# This converts real values to the wanted precision.
convert_floats(::Type{T}, g::GaussianSplitting{T}) where {T} = g

function convert_floats(::Type{T}, g::GaussianSplitting{S, N}) where {T, S, N}
    (; Ls, α, rcut, Ns) = g
    GaussianSplitting(convert.(T, Ls), convert(T, α), convert(T, rcut), Ns)
end

periods(g::GaussianSplitting) = g.Ls
cutoff_distance(g::GaussianSplitting) = g.rcut
fourier_grid_size(g::GaussianSplitting) = g.Ns

function _gaussian_splitting_params(Ls::NTuple{N, T}, β::Real, α::Real, rcut::Nothing, Ns::Nothing) where {N, T}
    rcut = β / α
    kmax = 2 * α * β
    Ns = map(Ls) do L
        kmax_to_gridsize(kmax, L)
    end
    T(α), T(rcut), Ns
end

function _gaussian_splitting_params(Ls::NTuple{N, T}, β::Real, α::Nothing, rcut::Real, Ns::Nothing) where {N, T}
    α = β / rcut
    kmax = 2 * α * β
    Ns = map(Ls) do L
        kmax_to_gridsize(kmax, L)
    end
    T(α), T(rcut), Ns
end

function _gaussian_splitting_params(Ls::NTuple{N, T}, β::Real, α::Nothing, rcut::Nothing, Ns::Dims{N}) where {N, T}
    kmax = maximum_wavenumber(Ns, Ls)
    α = kmax / (2 * β)
    rcut = β / α
    T(α), T(rcut), Ns
end

function _gaussian_splitting_params(Ls::NTuple{N, T}, β::Nothing, α::Real, rcut::Real, Ns::Dims{N}) where {N, T}
    T(α), T(rcut), Ns
end

function Base.show(io::IO, g::GaussianSplitting{T, N}) where {T, N}
    (; Ns, Ls, α, rcut) = g
    indent = get(io, :indent, 0)
    pre = ' '^indent
    rcut_L = rcut / minimum(Ls)
    kmax = maximum_wavenumber(Ns, Ls)
    β_shortrange = α * rcut
    β_longrange = kmax / (2 * α)
    print(io, "$(pre)GaussianSplitting{$T, $N} with:")
    print(io, "\n$(pre) - Domain period:               Ls = ", Ls)
    print(io, "\n$(pre) - Ewald splitting parameter:   α  = ", α)
    print(io, "\n$(pre) - Short-range cut-off:         r_cut = ", rcut, " (r_cut/L_min = ", rcut_L, ")")
    print(io, "\n$(pre) - Long-range resolution:       Ns = ", Ns, " (k_max = ", kmax, ")")
    print(io, "\n$(pre) - Short-range accuracy coeff.: β_shortrange = α * r_cut = ", β_shortrange)
    print(io, "\n$(pre) - Long-range accuracy coeff.:  β_longrange = k_max / 2α = ", β_longrange)
    nothing
end

function Base.summary(io::IO, g::GaussianSplitting)
    (; Ns) = g
    Ls = round.(g.Ls; sigdigits = 3)
    α = round(g.α; sigdigits = 3)
    rcut = round(g.rcut; sigdigits = 3)
    print(io, "GaussianSplitting(Ls ≈ $Ls, rcut ≈ $rcut, Ns = $Ns, α ≈ $α")
end

# Evaluate splitting kernel in Fourier space.
# Note that this may be called from a GPU kernel.
# It doesn't need to be very performant since it's only done once when creating a BiotSavartCache.
function splitting_kernel_fourier(g::GaussianSplitting)
    (; α,) = g
    β = -1 / (4 * α * α)
    @inline(k² -> exp(β * k²))
end

# Factor in ⟨ψ⟩ = C * ⟨ω⟩ to be applied when the mean vorticity ⟨ω⟩ is nonzero.
# See background_vorticity_correction! for details.
background_vorticity_correction_factor(g::GaussianSplitting) = 1 / (4 * g.α^2)

include("SIMDFunctions/SIMDFunctions.jl")
using .SIMDFunctions: SIMDFunctions

@inline exp_simd(x::SIMD.Vec) = SIMDFunctions.exp(x)
@inline erf_simd(x::SIMD.Vec) = SIMDFunctions.erf(x)
@inline erfc_simd(x::SIMD.Vec) = SIMDFunctions.erfc(x)

# Note: even without explicit SIMD, calling SIMD-friendly implementations can enable
# automatic SIMD and thus noticeably improve performance.
@inline exp_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.exp(x)
@inline erf_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.erf(x)
@inline erfc_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.erfc(x)

# On GPU we call the functions from Base or SpecialFunctions, since these are usually
# overridden in each GPU implementation (CUDA, ...) and therefore should be fast.
@inline exp_nosimd(::GPU, x::AbstractFloat) = exp(x)
@inline erf_nosimd(::GPU, x::AbstractFloat) = SpecialFunctions.erf(x)
@inline erfc_nosimd(::GPU, x::AbstractFloat) = SpecialFunctions.erfc(x)

@inline two_over_sqrt_pi(::SIMD.Vec{W, T}) where {W, T} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))

@inline function weights_shortrange_simd(g::GaussianSplitting, r)
    (; α,) = g
    αr = α * r
    a = erfc_simd(αr)
    b = two_over_sqrt_pi(r) * αr * exp_simd(-(αr * αr))
    a, b
end

@inline function weights_shortrange_nosimd(backend::KA.Backend, g::GaussianSplitting, r)
    (; α,) = g
    αr = α * r
    a = erfc_nosimd(backend, αr)
    b = two_over_sqrt_pi(r) * αr * exp_nosimd(backend, -(αr * αr))
    a, b
end

# Same as above but replacing erfc -> erf.
# @inline function weights_longrange_simd(g::GaussianSplitting, r)
#     (; α,) = g
#     αr = α * r
#     a = erf_simd(αr)
#     b = two_over_sqrt_pi(r) * αr * exp_simd(-(αr * αr))
#     a, b
# end

@inline function weights_longrange_nosimd(backend::KA.Backend, g::GaussianSplitting, r)
    (; α,) = g
    αr = α * r
    a = erf_nosimd(backend, αr)
    b = two_over_sqrt_pi(r) * αr * exp_nosimd(backend, -(αr * αr))
    a, b
end
