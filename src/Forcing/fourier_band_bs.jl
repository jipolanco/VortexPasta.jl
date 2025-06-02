using ..BiotSavart: BiotSavart, BiotSavartCache
using StaticArrays: SVector, SMatrix
using LinearAlgebra: LinearAlgebra, ⋅, ×

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
\begin{align*}
\frac{δE(\bm{k})}{δ\bm{s}}
&= \frac{Γ}{2V} i \bm{B} \bm{ψ}(\bm{k}) + \text{c.c.}
\\
&= \frac{Γ}{2V} [\bm{s}' × \bm{v}(\bm{k})] e^{i \bm{k} ⋅ \bm{s}} + \text{c.c.}
\\
&= ℜ \left\{ \frac{Γ}{V} [\bm{s}' × \bm{v}(\bm{k})] e^{i \bm{k} ⋅ \bm{s}} \right\}
\end{align*}
```

Finally, the idea is to advect the filaments with a velocity which is parallel to this result for each forced ``\bm{k}``.
One can write this velocity as

```math
\bm{v}_{\text{f}}(\bm{k}, \bm{s}) = α \, ℜ \left\{ [\bm{s}' × \bm{v}(\bm{k})] e^{i \bm{k} ⋅ \bm{s}} \right\} = α \, \bm{v}_0(\bm{k}, \bm{s})
```

where ``α`` is a non-dimensional parameter setting the forcing amplitude.

## Estimating the energy injection rate

One can also try to estimate an energy injection rate at wavevector ``\bm{k}`` associated to this velocity:

```math
\frac{\text{d}E(\bm{k})}{\text{d}t}
= ∮ \frac{δE(\bm{k})}{δ\bm{s}} ⋅ \frac{\text{d}\bm{s}}{\text{d}t} \, \mathrm{d}ξ
= α \frac{Γ}{V} ∮ |\bm{v}_0|^2 \, \mathrm{d}ξ
```

In general, this estimate may be quite inaccurate since the forcing can also affect the
energy at wavevectors other than ``\bm{k}``. But still, this estimate is the one used when
``ε_{\text{target}}`` is given, and may allow to obtain a roughly constant energy injection
rate (even if it's generally different than the "target" one).
"""
struct FourierBandForcingBS{T <: AbstractFloat, N} <: AbstractForcing
    α    :: T
    ε_target :: T  # target energy injection rate [L²T⁻³]
    kmin :: T
    kmax :: T
    qs   :: Vector{NTuple{N, Int}}  # list of normalised wavevectors (can be empty)
end

function FourierBandForcingBS(; α::Real = 0, ε_target = 0, kmin = nothing, kmax = nothing, qs = nothing)
    (α == 0) + (ε_target == 0) == 1 || throw(ArgumentError("one should pass either α or ε_target, but not both"))
    _FourierBandForcingBS(α, ε_target, kmin, kmax, qs)
end

function _FourierBandForcingBS(α, ε_target, kmin::Real, kmax::Real, qs::Nothing)
    scalars = promote(α, ε_target, kmin, kmax)
    qs = Tuple{}[]  # empty vector of empty tuples
    FourierBandForcingBS(scalars..., qs)
end

# This variant allows setting specific q⃗ normalised wavevectors (currently not documented!).
function _FourierBandForcingBS(α, ε_target, kmin::Nothing, kmax::Nothing, qs::AbstractVector{NTuple{N, Int}}) where {N}
    kmin = 0
    kmax = 0
    scalars = promote(α, ε_target, kmin, kmax)
    FourierBandForcingBS(scalars..., qs)
end

function Base.show(io::IO, f::FourierBandForcingBS{T}) where {T}
    (; α, ε_target, kmin, kmax, qs,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "FourierBandForcingBS{$T} with:")
    if α != 0
        print(io, "\n$(prefix)├─ Magnitude: α = ", α)
    elseif ε_target != 0
        print(io, "\n$(prefix)├─ Target energy injection rate: ε_target = ", ε_target)
    end
    if isempty(qs)
        print(io, "\n$(prefix)└─ Fourier band: |k⃗| ∈ [$kmin, $kmax]")
    else
        print(io, "\n$(prefix)└─ Normalised Fourier wave vectors: |q⃗| = ", qs)
    end
end

# Here vs_grid is the superfluid velocity in Fourier space, optionally on a device (GPU).
function init_cache(f::FourierBandForcingBS, cache_bs::BiotSavartCache)
    (; kmin, kmax, qs,) = f
    (; params,) = cache_bs
    vs_grid = BiotSavart.get_longrange_field_fourier(cache_bs).field :: NTuple
    backend = KA.get_backend(vs_grid[1])  # CPU, CUDABackend, ROCBackend, ...
    if isempty(qs)
        # Activate all wavevectors within a Fourier band [kmin, kmax]
        v_h = FourierBandVectorField(undef, params.Ls; kmin, kmax)  # on the host (CPU)
    else
        # Activate selected wavenumbers
        v_h = FourierBandVectorField(undef, params.Ls)  # on the host (CPU)
        for q⃗ in qs
            SyntheticFields.add_normalised_wavevector!(v_h, q⃗)
        end
    end
    v_d = adapt(backend, v_h)  # on the "device" (GPU or CPU)
    prefactor = params.Γ / prod(params.Ls)
    (; v_d, v_h, qs = v_h.qs, prefactor,)
end

# After calling this function, cache.v_h contains the unfiltered Biot-Savart velocity in Fourier space v̂(k⃗).
function update_cache!(cache, f::FourierBandForcingBS, cache_bs::BiotSavartCache)
    (; v_d, v_h,) = cache

    vs_grid, ks_grid, σ_gaussian = let data = BiotSavart.get_longrange_field_fourier(cache_bs)
        local (; state, field, wavenumbers,) = data
        @assert state.quantity == :velocity
        field, wavenumbers, state.smoothing_scale
    end

    # (1) Compute unfiltered velocity within Fourier band
    σ²_over_two = σ_gaussian^2 / 2
    @inline function op_streamfunction(_, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        φ = exp(k² * σ²_over_two)
        φ * vs_filtered
    end
    SyntheticFields.from_fourier_grid!(op_streamfunction, v_d, vs_grid, ks_grid)

    # (2) Copy results to CPU if needed (avoided if the "device" is the CPU).
    if v_d !== v_h
        copyto!(v_h, v_d)
    end

    nothing
end

# Returns the estimated energy injection rate.
function apply!(
        forcing::FourierBandForcingBS, cache, vs::AbstractVector, f::AbstractFilament;
        scheduler = SerialScheduler(),
    )
    (; v_h,) = cache
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    ts = Filaments.knots(f)
    ε_total = tmapreduce(+, eachindex(vs); scheduler) do i
        @inline
        s⃗ = f[i]
        ds⃗_dt = f[i, Derivative(1)]
        ds_dt = sqrt(sum(abs2, ds⃗_dt))  # vector norm
        s⃗′ = ds⃗_dt ./ ds_dt  # unit tangent
        (; qs, cs, Δks,) = v_h
        vf = zero(V)  # forcing velocity (excluding α prefactor)
        for n in eachindex(qs, cs)  # iterate over active wavevectors k⃗
            k⃗ = SVector(qs[n]) .* Δks
            v̂ = cs[n]
            vf_k = s⃗′ × real(v̂ * cis(k⃗ ⋅ s⃗))  # forcing velocity for this wavevector k⃗
            vf = vf + vf_k
        end
        if forcing.α != 0
            vs[i] = vs[i] + forcing.α * vf  # apply forcing
        elseif forcing.ε_target != 0
            vs[i] = vf  # just overwrite the input
        end
        dt = (ts[i + 1] - ts[i - 1]) / 2
        dξ = dt * ds_dt  # estimated local segment length
        sum(abs2, vf) * dξ  # estimated energy injection rate around local node (without Γ/V prefactor)
    end
    # The factor 2 accounts for Hermitian symmetry: if we inject energy into the velocity at
    # a wavevector +k⃗, then we're also injecting the same energy at v(-k⃗).
    2 * ε_total * cache.prefactor  # include Γ/V prefactor
end
