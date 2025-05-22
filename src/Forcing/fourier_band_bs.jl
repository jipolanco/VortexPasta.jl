using ..BiotSavart: BiotSavart, BiotSavartCache
using StaticArrays: SVector, SMatrix
using LinearAlgebra: LinearAlgebra, ⋅

@doc raw"""
    FourierBandForcingBS <: AbstractForcing
    FourierBandForcingBS(; kmin, kmax, α,)

Forcing based on Biot–Savart energetics within a given range of wavenumbers.

This type of forcing does _not_ rely on an imposed normal fluid velocity. Instead, it starts from the
functional derivative of the energy at a given wavevector ``\bm{k}`` with respect to the
vortex positions (according to the Biot–Savart law), and applies a velocity that ensures
energy to increase (if ``α > 0``) at those wavevectors.
"""
struct FourierBandForcingBS{T <: AbstractFloat} <: AbstractForcing
    α    :: T
    kmin :: T
    kmax :: T
end

function FourierBandForcingBS(; α::Real, kmin::Real, kmax::Real)
    FourierBandForcingBS(promote(α, kmin, kmax)...)
end

function Base.show(io::IO, f::FourierBandForcingBS{T}) where {T}
    (; α, kmin, kmax,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "FourierBandForcingBS{$T} with:")
    print(io, "\n$(prefix)├─ Magnitude: α = ", α)
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
    (; ψ_d, ψ_h,)
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
    u = similar(ψ̂)
    for i in eachindex(u)
        u[i] = e * (k⃗[i] * ψs′ - ψ̂[i] * ks′)
    end
    SVector(u)
end

function apply!(
        forcing::FourierBandForcingBS, cache, vs::AbstractVector, f::AbstractFilament;
        scheduler = SerialScheduler(),
    )
    (; α,) = forcing
    (; ψ_h,) = cache
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    tforeach(eachindex(vs); scheduler) do i
        @inline
        s⃗ = f[i]
        s⃗′ = f[i, UnitTangent()]
        vf = zero(V)  # forcing velocity (excluding α prefactor)
        (; qs, cs, Δks,) = ψ_h
        for n in eachindex(qs, cs)  # iterate over active wavevectors k⃗
            k⃗ = SVector(qs[n]) .* Δks
            ψ̂ = cs[n]
            u = _apply_forcing_matrix(k⃗, s⃗, s⃗′, ψ̂)
            vf = vf - imag(u)
        end
        vs[i] = vs[i] + α * vf  # apply forcing
    end
    nothing
end
