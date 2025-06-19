@doc raw"""
    SmallScaleDissipationBS <: AbstractDissipation
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
struct SmallScaleDissipationBS{T <: AbstractFloat} <: AbstractDissipation
    kdiss :: T
    α     :: T
    ε_target :: T
end

function SmallScaleDissipationBS(; kdiss::Real, α::Real = 0, ε_target::Real = 0)
    (α == 0) + (ε_target == 0) == 1 || throw(ArgumentError("one should pass either α or ε_target, but not both"))
    kdiss ≥ 0 || throw(ArgumentError("dissipative wavenumber `kdiss` should be non-negative"))
    SmallScaleDissipationBS(promote(kdiss, α, ε_target)...)
end

function Base.show(io::IO, f::SmallScaleDissipationBS{T}) where {T}
    (; α, ε_target, kdiss,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "SmallScaleDissipationBS{$T} with:")
    if α != 0
        print(io, "\n$(prefix)├─ Magnitude: α = ", α)
    elseif ε_target != 0
        print(io, "\n$(prefix)├─ Target energy dissipation rate: ε_target = ", ε_target)
    end
    print(io, "\n$(prefix)└─ Fourier wavenumbers: |k⃗| > $kdiss")
    nothing
end

function init_cache(dissipation::SmallScaleDissipationBS, cache_bs::BiotSavartCache)
    (; kdiss,) = dissipation
    (; params,) = cache_bs
    # Check that the dissipation wavenumber is inside the Fourier grid in all 3 directions
    (; wavenumbers,) = BiotSavart.get_longrange_field_fourier(cache_bs)
    for ks in wavenumbers  # iterate over the 3 dimensions
        if ks[end] > 0  # r2c FFT, usually the first dimension
            kmax = ks[end - 1]  # discard the very last wavenumber
        else  # c2c FFT
            kmax = ks[(end + 1) ÷ 2]
            @assert kmax > 0
        end
        if kdiss ≥ kmax
            error(lazy"SmallScaleDissipationBS: dissipative wavenumber (kdiss = $kdiss) is too large compared to size of Fourier grid (kmax = $kmax)")
        end
    end
    prefactor = params.Γ / prod(params.Ls)
    (; longrange = cache_bs.longrange, prefactor,)
end

function _update_cache!(cache, dissipation::SmallScaleDissipationBS)
    # Interpolate Fourier-truncated Biot-Savart velocities at filament locations (results
    # are written somewhere in cache.longrange).
    (; kdiss,) = dissipation
    (; wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache.longrange)
    @assert state.quantity == :velocity && state.smoothing_scale == 0
    # Interpolate Fourier velocity, but first discarding Fourier coefficients k > kdiss.
    kdiss² = kdiss * kdiss
    BiotSavart.interpolate_to_physical!(cache.longrange) do û::SVector{3}, I::CartesianIndex{3}
        @inline
        @inbounds k⃗ = getindex.(wavenumbers, Tuple(I))
        k² = sum(abs2, k⃗)
        ifelse(k² > kdiss², zero(û), û)
    end
    nothing
end

function apply!(
        dissipation::SmallScaleDissipationBS{T}, cache,
        vL_all::AbstractVector, vdiss_all::AbstractVector, vs_all::AbstractVector, fs::AbstractVector{<:AbstractFilament};
        scheduler = SerialScheduler(),
    ) where {T <: AbstractFloat}
    # 1. Interpolate vs^< (Fourier-truncated BS velocity) onto cache.longrange
    _update_cache!(cache, dissipation)

    # 2. Copy results onto vdiss_all
    BiotSavart.copy_long_range_output!(vdiss_all, cache.longrange)

    # 3. Iterate over filaments and filament points to obtain the actual dissipation velocity.
    ε_total = zero(T)
    @inbounds for n in eachindex(fs)
        f = fs[n]
        vs = vs_all[n]  # Biot-Savart velocity of filament
        vdiss = vdiss_all[n]  # for now this contains vs^< on the filament
        ts = Filaments.knots(f)
        ε_filament = tmapreduce(+, eachindex(vs); scheduler) do i
            @inline
            @inbounds begin
                ds⃗_dt = f[i, Derivative(1)]
                ds_dt = sqrt(sum(abs2, ds⃗_dt))  # vector norm
                s⃗′ = ds⃗_dt ./ ds_dt  # unit tangent
                v⃗_hi = vs[i] - vdiss[i]  # high-pass filtered velocity
                v⃗_diss = v⃗_hi × s⃗′  # dissipation velocity (without α coefficient) -- this already includes the negative sign for dissipation: v⃗ × s⃗′ = -s⃗′ × v⃗
                if dissipation.α != 0
                    vdiss[i] = dissipation.α * v⃗_diss
                else
                    vdiss[i] = v⃗_diss  # we will rescale vdiss later to get the wanted dissipation rate
                end
                dt = (ts[i + 1] - ts[i - 1]) / 2
                dξ = dt * ds_dt  # estimated local segment length
                sum(abs2, v⃗_diss) * dξ  # estimated energy dissipation rate around local node (without Γ/V prefactor) -- only valid in ε_target mode!
            end
        end
        ε_total += ε_filament
    end
    ε_total *= cache.prefactor

    # 4. Rescale velocities if needed
    if dissipation.ε_target != 0 && ε_total != 0
        @assert dissipation.α == 0
        # TODO: check that abs(ε_total) > some minimal tolerance? (to avoid division by "almost" zero)
        α = dissipation.ε_target / ε_total
        @. vdiss_all *= α
    end

    # 5. Add dissipation term to vortex velocity vL
    @. vL_all = vL_all + vdiss_all

    nothing
end
