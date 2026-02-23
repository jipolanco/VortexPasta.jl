@doc raw"""
    KaiserBesselSplitting{T <: AbstractFloat} <: AbstractEwaldSplitting
    KaiserBesselSplitting(; Ls, β, rcut, Ns)

Kaiser--Bessel splitting kernel for Ewald summation.

This kernel presents near-optimal localisation in Fourier space, meaning that, for a given
cut-off distance `rcut` in physical space, one can use a smaller Fourier grid size `Ns` to
achieve the same accuracy as [`GaussianSplitting`](@ref).

Note that, unlike `GaussianSplitting`, this kernel has zero truncation error in physical space,
since the Kaiser--Bessel kernel has compact support (`|r| < rcut`).

# Parameters

- `Ls::NTuple{3, T}`: domain period in each Cartesian direction (**mandatory**),
- `β::T`: nondimensional shape parameter,
- `rcut::T`: cut-off distance for short-range computations,
- `Ns::NTuple{3, Int}`: size of FFT grid for long-range computations.

Note that one does _not_ need to pass all of these parameters (see below).

# Construction

There are several ways of constructing a `KaiserBesselSplitting`.
Note that the period `Ls` is always required.

The **recommended way** is to pass the nondimensional accuracy coefficient `β` and _one of_
`rcut` or `Ns`. For example:

```jldoctest
julia> Ls = (2π, 2π, 2π);
julia> β = 14;  # gives roughly 6-digit accuracy (to be checked...)
julia> Ns = (128, 128, 128);
julia> splitting = KaiserBesselSplitting(; Ls, β, Ns)
KaiserBesselSplitting{Float64, 3} with:
 - Domain period:          Ls = (6.283185307179586, 6.283185307179586, 6.283185307179586)
 - Shape parameter:        β  = 14.0
 - Short-range cut-off:    r_cut = 0.2222222222222222 (r_cut/L_min = 0.035367765131532294)
 - Long-range resolution:  Ns = (128, 128, 128) (k_max = 63.0)
```

# Extended help

## Kernel definitions

This kernel splits the 3D Green's function ``G(\bm{r}) = 1 / (4πr)`` into the near- and far-field
contributions:

```math
G(\bm{r}) = G^{\text{(n)}}(\bm{r}) + G^{\text{(f)}}(\bm{r}) =
\frac{1 - F(r)}{4πr} + \frac{F(r)}{4πr}
```

where

```math
F(r) = 2 \int_0^r f(u) \, \mathrm{d}u
```

and ``f(r)`` is the Kaiser--Bessel kernel:

```math
f(r) = \frac{β}{2 r_{\text{c}} \sinh(β)} I_0 \! \left(β \sqrt{1 - \frac{r^2}{r_{\text{c}}^2}}\right)
\quad \text{for } |r| ≤ r_{\text{c}},
```

normalised such that ``F(r_{\text{c}}) = 1``.
Here ``r_{\text{c}}`` represents the kernel support in physical space while ``β`` is a
nondimensional shape parameter.

As a result, the Biot--Savart kernel $\bm{\nabla}G(\bm{r}) = -\bm{r} / (4πr^3)$ is split as
$\bm{\nabla}G(\bm{r}) = \bm{\nabla}G^{\text{(n)}}(\bm{r}) + \bm{\nabla}G^{\text{(f)}}(\bm{r})$ with:

```math
\begin{align*}
    \bm{\nabla}G^{\text{(n)}}(\bm{r})
    &= -\frac{\bm{r}}{4\pi r^3} \left[ 1 - F(r) + 2r f(r) \right]
    \\
    \bm{\nabla}G^{\text{(f)}}(\bm{r})
    &= -\frac{\bm{r}}{4\pi r^3} \left[ F(r) - 2r f(r) \right]
\end{align*}
```
"""
struct KaiserBesselSplitting{T <: AbstractFloat, N} <: AbstractEwaldSplitting
    Ls::NTuple{N, T}
    β::T
    rcut::T
    Ns::Dims{N}
end

function KaiserBesselSplitting(; Ls::NTuple{N, T}, β = nothing, rcut = nothing, Ns = nothing) where {N, T}
    KaiserBesselSplitting(Ls, _kb_splitting_params(Ls, β, rcut, Ns)...)
end

# This converts real values to the wanted precision.
convert_floats(::Type{T}, g::KaiserBesselSplitting{T}) where {T} = g

function convert_floats(::Type{T}, g::KaiserBesselSplitting{S, N}) where {T, S, N}
    (; Ls, β, rcut, Ns) = g
    KaiserBesselSplitting(convert.(T, Ls), convert(T, β), convert(T, rcut), Ns)
end

periods(g::GaussianSplitting) = g.Ls
cutoff_distance(g::GaussianSplitting) = g.rcut
fourier_grid_size(g::GaussianSplitting) = g.Ns

# β and rcut given
function _kb_splitting_params(Ls::NTuple{N, T}, β::Real, rcut::Real, Ns::Nothing) where {N, T}
    kmax = β / rcut
    Ns = map(Ls) do L
        m = floor(Int, kmax * L / T(2π))
        2m + 1
    end
    T(β), T(rcut), Ns
end

# rcut and Ns given
function _kb_splitting_params(Ls::NTuple{N, T}, β::Nothing, rcut::Real, Ns::Dims{N}) where {N, T}
    kmax = maximum_wavenumber(Ns, Ls)
    β = rcut * kmax
    T(β), T(rcut), Ns
end

# β and Ns given
function _kb_splitting_params(Ls::NTuple{N, T}, β::Real, rcut::Nothing, Ns::Dims{N}) where {N, T}
    kmax = maximum_wavenumber(Ns, Ls)
    rcut = β / kmax
    T(β), T(rcut), Ns
end

function Base.show(io::IO, g::KaiserBesselSplitting{T, N}) where {T, N}
    (; Ns, Ls, β, rcut) = g
    indent = get(io, :indent, 0)
    pre = ' '^indent
    rcut_L = rcut / minimum(Ls)
    kmax = maximum_wavenumber(Ns, Ls)
    print(io, "$(pre)KaiserBesselSplitting{$T, $N} with:")
    print(io, "\n$(pre) - Domain period:          Ls = ", Ls)
    print(io, "\n$(pre) - Shape parameter:        β  = ", β)
    print(io, "\n$(pre) - Short-range cut-off:    r_cut = ", rcut, " (r_cut/L_min = ", rcut_L, ")")
    print(io, "\n$(pre) - Long-range resolution:  Ns = ", Ns, " (k_max = ", kmax, ")")
    nothing
end
