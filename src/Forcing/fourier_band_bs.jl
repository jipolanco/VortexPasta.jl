using ..BiotSavart: BiotSavart, BiotSavartCache
using StaticArrays: SVector, SMatrix
using LinearAlgebra: LinearAlgebra, ⋅

@doc raw"""
    FourierBandForcingBS <: AbstractForcing
    FourierBandForcingBS(; kmin, kmax, α, ε_target)

Forcing based on Biot–Savart energetics within a given range of wavenumbers.

This type of forcing does _not_ rely on an imposed normal fluid velocity. Instead, it starts from the
functional derivative of the energy at a given wavevector ``\bm{k}`` with respect to the
vortex positions (according to the Biot–Savart law), and applies a velocity that ensures
energy to increase (if ``α > 0`` or ``ε_{\text{target}} > 0``) at those wavevectors.

One should pass _either_ `α` or `ε_target` but never both. They should be positive for
energy injection (negative values lead to energy dissipation):

- `α` (`\alpha`) is a non-dimensional coefficient which directly sets the amplitude of the forcing velocity;

- `ε_target` (`\varepsilon_target`) has the units of an energy injection rate. In this case,
  the amplitude ``α`` will be adjusted over time in order to keep a roughly constant energy
  injection rate (which in general will _not_ be equal to `ε_target`, see remarks below).

# Extended help

## Forcing definition

This forcing attempts to increase the kinetic energy at selected wavenumbers.
Its definition starts from the expression for the kinetic energy at wavenumber ``\bm{k}``:

```math
E(\bm{k}) = \frac{1}{2} |\bm{v}(\bm{k})|^2 = \frac{1}{2k^2} |\bm{ω}(\bm{k})|^2
```

The idea is to translate the vortex positions by ``\bm{s}(ξ) → \bm{s}(ξ) + δ\bm{s}(ξ)`` so
that the energy ``E(\bm{k})`` increases. To determine such a displacement, one can look at
the functional derivative of ``E(\bm{k})`` with respect to the positions ``\bm{s}``:

```math
\begin{align*}
\frac{δE(\bm{k})}{δ\bm{s}}
&= \frac{1}{2k^2} \bm{ω}(\bm{k}) ⋅ \frac{δ\bm{ω}^*(\bm{k})}{δ\bm{s}} + \text{c.c.},
\\
&= \frac{1}{2} \bm{ψ}(\bm{k}) ⋅ \frac{δ\bm{ω}^*(\bm{k})}{δ\bm{s}} + \text{c.c.},
\end{align*}
```

where ``()^*`` denotes a complex conjugate and "c.c." means the complex conjugate of the first term.
Here we have used ``\bm{ω}(\bm{k}) = k^2 \bm{ψ}(\bm{k})``, which corresponds in physical
space to ``\bm{ω} = -∇² \bm{ψ}``.

The functional derivative of vorticity in Fourier space is:

```math
\begin{align*}
\frac{δ\bm{ω}^*(\bm{k})}{δ\bm{s}}
&= \frac{Γ}{V} \frac{δ}{δ\bm{s}} ∮ e^{+i \bm{k} ⋅ \bm{s}(ξ)} \bm{s}(ξ)' \, \text{d}ξ
\\
&= \frac{Γ}{V} \frac{δ}{δ\bm{s}} ∮ \bm{ℒ}[\bm{s}(ξ), \bm{s}(ξ)'] \, \text{d}ξ
\\
&= \frac{Γ}{V} \left[ \frac{∂\bm{ℒ}}{∂\bm{s}} - \frac{\text{d}}{\text{d}ξ} \frac{∂\bm{ℒ}}{∂\bm{s}'} \right]
\\
&= \frac{i Γ}{V} e^{+i \bm{k} ⋅ \bm{s}} \left[ \bm{k} ⊗ \bm{s}' - (\bm{k} ⋅ \bm{s}') I \right] = \frac{i Γ}{V} \bm{B}(\bm{k}, \bm{s}, \bm{s}')
\end{align*}
```

where ``\bm{B}`` is a ``3×3`` complex matrix which has dimensions of an inverse length (``L^{-1}``).

The functional derivative of ``E(\bm{k})`` can now be written as:

```math
\frac{δE(\bm{k})}{δ\bm{s}}
= \frac{iΓ}{2V} \bm{B} \bm{ψ}(\bm{k}) + \text{c.c.}
= ℜ \left\{ \frac{iΓ}{V} \bm{B} \bm{ψ}(\bm{k}) \right\}
= -ℑ \left\{ \frac{Γ}{V} \bm{B} \bm{ψ}(\bm{k}) \right\}
```

Finally, the idea is to advect the filaments with a velocity which is parallel to this result for each forced ``\bm{k}``.
One can write this velocity as

```math
\bm{v}_{\text{f}}(\bm{k}, \bm{s}) = α \, \Re\left\{ i \bm{B} \bm{ψ}(\bm{k}) \right\} = α \, \bm{v}_0(\bm{k}, \bm{s})
```

where ``α`` is a non-dimensional parameter setting the forcing amplitude.

## Estimating the energy injection rate

One can also try to estimate an energy injection rate at wavevector ``\bm{k}`` associated to this velocity:

```math
\frac{\text{d}E(\bm{k})}{\text{d}t}
= ∮ \frac{δE(\bm{k})}{δ\bm{s}} \, \frac{\text{d}\bm{s}}{\text{d}t} \, \mathrm{d}ξ
= α \frac{Γ}{V} ∮ |\bm{v}_0|^2 \, \mathrm{d}ξ
```

Unfortunately, this estimate is generally quite inaccurate since the forcing can also affect
the energy at wavevectors other than ``\bm{k}``. Conversely, ``E(\bm{k})`` may also be
affected by the forcing velocity associated to other forced wavevectors. But still, this
estimate is the one used when ``ε_{\text{target}}`` is given, and may allow to obtain a
roughly constant energy injection rate (even if it's generally different than the "target" one).
"""
struct FourierBandForcingBS{T <: AbstractFloat} <: AbstractForcing
    α    :: T
    ε_target :: T  # target energy injection rate [L²T⁻³]
    kmin :: T
    kmax :: T
end

function FourierBandForcingBS(; α::Real = 0, ε_target = 0, kmin::Real, kmax::Real)
    α * ε_target == 0 && α + ε_target != 0 || throw(ArgumentError("one should pass either α or ε_target, but not both"))
    FourierBandForcingBS(promote(α, ε_target, kmin, kmax)...)
end

function Base.show(io::IO, f::FourierBandForcingBS{T}) where {T}
    (; α, ε_target, kmin, kmax,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "FourierBandForcingBS{$T} with:")
    if α != 0
        print(io, "\n$(prefix)├─ Magnitude: α = ", α)
    elseif ε_target != 0
        print(io, "\n$(prefix)├─ Target energy injection rate: ε_target = ", ε_target)
    end
    print(io, "\n$(prefix)└─ Fourier band: |k⃗| ∈ [$kmin, $kmax]")
end

# Here vs_grid is the superfluid velocity in Fourier space, optionally on a device (GPU).
function init_cache(f::FourierBandForcingBS, cache_bs::BiotSavartCache)
    (; kmin, kmax,) = f
    (; params,) = cache_bs
    vs_grid = BiotSavart.get_longrange_field_fourier(cache_bs).field :: NTuple
    backend = KA.get_backend(vs_grid[1])  # CPU, CUDABackend, ROCBackend, ...
    ψ_h = FourierBandVectorField(undef, params.Ls; kmin, kmax)  # on the host (CPU)
    ψ_d = adapt(backend, ψ_h)  # on the "device" (GPU or CPU)
    prefactor = params.Γ / prod(params.Ls)
    (; ψ_d, ψ_h, prefactor,)
end

# After calling this function, cache.ψ_h contains ψ̂(k⃗) where ψ is the streamfunction.
function update_cache!(cache, f::FourierBandForcingBS, cache_bs::BiotSavartCache)
    (; ψ_d, ψ_h,) = cache

    vs_grid, ks_grid = let data = BiotSavart.get_longrange_field_fourier(cache_bs)
        local (; state, field, wavenumbers,) = data
        @assert state.quantity == :velocity
        @assert state.smoothed == true
        field, wavenumbers
    end
    α_ewald = cache_bs.params.α

    # (1) Compute unfiltered streamfunction from filtered velocity within Fourier band
    inv_four_α² = 1 / (4 * α_ewald * α_ewald)
    @inline function op_streamfunction(_, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        φ = @fastmath exp(k² * inv_four_α²)
        vs = φ * vs_filtered
        im * (k⃗ × vs) ./ k²
    end
    SyntheticFields.from_fourier_grid!(op_streamfunction, ψ_d, vs_grid, ks_grid)

    # (2) Copy results to CPU if needed (avoided if the "device" is the CPU).
    if ψ_d !== ψ_h
        copyto!(ψ_h, ψ_d)
    end

    nothing
end

@inline function _apply_forcing_matrix(k⃗::SVector{N}, s⃗::SVector{N}, s⃗′::SVector{N}, ψ̂::SVector{N}) where {N}
    ks′ = k⃗ ⋅ s⃗′
    ψs′ = s⃗′ ⋅ ψ̂    # note: ψ̂ must be on the right because it's complex (otherwise it will be conjugated)
    e = cis(k⃗ ⋅ s⃗)  # = exp(im * (k⃗ ⋅ s⃗))
    e * (k⃗ * ψs′ - ψ̂ * ks′)
end

# Returns the estimated energy injection rate (a very rough estimate likely to be wrong, but
# may still be useful to obtain a steady state).
function apply!(
        forcing::FourierBandForcingBS, cache, vs::AbstractVector, f::AbstractFilament;
        scheduler = SerialScheduler(),
    )
    (; ψ_h,) = cache
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    ts = Filaments.knots(f)
    ε_total = tmapreduce(+, eachindex(vs); scheduler) do i
        @inline
        s⃗ = f[i]
        ds⃗_dt = f[i, Derivative(1)]
        ds_dt = sqrt(sum(abs2, ds⃗_dt))  # vector norm
        s⃗′ = ds⃗_dt ./ ds_dt  # unit tangent
        (; qs, cs, Δks,) = ψ_h
        dt = (ts[i + 1] - ts[i - 1]) / 2
        dξ = dt * ds_dt  # estimated local segment length
        ε = zero(eltype(s⃗))
        vf = zero(V)  # forcing velocity (excluding α prefactor)
        for n in eachindex(qs, cs)  # iterate over active wavevectors k⃗
            k⃗ = SVector(qs[n]) .* Δks
            ψ̂ = cs[n]
            u = _apply_forcing_matrix(k⃗, s⃗, s⃗′, ψ̂)
            vf_k = -imag(u)  # forcing velocity for this wavevector k⃗
            ε += sum(abs2, vf_k)  # energy injection estimate (without prefactor Γ dξ / V)
            vf = vf + vf_k
        end
        if forcing.α != 0
            vs[i] = vs[i] + forcing.α * vf  # apply forcing
        elseif forcing.ε_target != 0
            vs[i] = vf  # just overwrite the input
        end
        ε * dξ  # estimated energy injection rate around local node
    end
    ε_total * cache.prefactor  # include Γ/V prefactor
end
