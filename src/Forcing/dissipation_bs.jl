@doc raw"""
    DissipationBS <: AbstractDissipation
    DissipationBS(; α, ε_target)

Dissipation based on Biot–Savart energetics.

Note that this dissipation method acts at **all scales**.
To dissipate at small scales only, use [`SmallScaleDissipationBS`](@ref).

To set the dissipation amplitude, one should pass _either_ `α` or `ε_target` but never both.
They should be positive for energy dissipation (negative values lead to energy injection,
which may lead to instabilities!):

- `α` (`\alpha`) is a non-dimensional coefficient which directly sets the amplitude of the
  dissipation term;

- `ε_target` (`\varepsilon_target`) has the units of an energy dissipation rate. In this case,
  the amplitude ``α`` will be adjusted over time in order to keep a constant energy
  dissipation rate.

# Extended help

## Dissipation term definition

The idea is to apply a "dissipation velocity" ``\bm{v}_{\text{diss}}`` to the vortices which
will extract kinetic energy from the system. Such dissipation term can be written as:

```math
\bm{v}_{\text{diss}}(\bm{s}) = -α \bm{s}' × \bm{v}(\bm{s})
```

where ``\bm{v}(\bm{s})`` is the superfluid velocity (from Biot–Savart's law).

## Energy dissipation rate

The dissipation rate associated to this term is:

```math
ε_{\text{diss}} = -α \frac{Γ}{V} ∮ |\bm{s}' × \bm{v}|^2 \, \mathrm{d}ξ
```
"""
struct DissipationBS{T <: AbstractFloat} <: AbstractDissipation
    α     :: T
    ε_target :: T
end

function DissipationBS(; α::Real = 0, ε_target::Real = 0)
    (α == 0) + (ε_target == 0) == 1 || throw(ArgumentError("one should pass either α or ε_target, but not both"))
    DissipationBS(promote(α, ε_target)...)
end

function Base.show(io::IO, f::DissipationBS{T}) where {T}
    (; α, ε_target,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "DissipationBS{$T} with:")
    if α != 0
        print(io, "\n$(prefix)└─ Magnitude: α = ", α)
    elseif ε_target != 0
        print(io, "\n$(prefix)└─ Target energy dissipation rate: ε_target = ", ε_target)
    end
    nothing
end

function init_cache(dissipation::DissipationBS, cache_bs::BiotSavartCache)
    (; params,) = cache_bs
    prefactor = params.Γ / prod(params.Ls)
    (; prefactor,)
end

function apply!(
        dissipation::DissipationBS{T}, cache,
        vL_all::AbstractVector, vdiss_all::AbstractVector, vs_all::AbstractVector, fs::AbstractVector{<:AbstractFilament};
        scheduler = SerialScheduler(),
    ) where {T <: AbstractFloat}
    # Iterate over filaments and filament points to obtain the actual dissipation velocity.
    ε_total = zero(T)
    @inbounds for n in eachindex(fs)
        f = fs[n]
        vs = vs_all[n]        # Biot-Savart velocity of filament
        vdiss = vdiss_all[n]  # dissipation velocity will be written here
        ts = Filaments.knots(f)
        ε_filament = tmapreduce(+, eachindex(vs); scheduler) do i
            @inline
            @inbounds begin
                ds⃗_dt = f[i, Derivative(1)]
                ds_dt = sqrt(sum(abs2, ds⃗_dt))  # vector norm
                s⃗′ = ds⃗_dt ./ ds_dt  # unit tangent
                v⃗ₛ = vs[i]
                v⃗_diss = v⃗ₛ × s⃗′  # dissipation velocity (without α coefficient) -- this already includes the negative sign for dissipation: v⃗ × s⃗′ = -s⃗′ × v⃗
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

    # Rescale velocities if needed
    if dissipation.ε_target != 0 && ε_total != 0
        @assert dissipation.α == 0
        # TODO: check that abs(ε_total) > some minimal tolerance? (to avoid division by "almost" zero)
        α = dissipation.ε_target / ε_total
        @. vdiss_all *= α
    end

    # Add dissipation term to vortex velocity vL
    @. vL_all = vL_all + vdiss_all

    nothing
end
