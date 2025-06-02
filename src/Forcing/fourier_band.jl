@doc raw"""
    FourierBandForcing <: AbstractForcing
    FourierBandForcing(vn::FourierBandVectorField; α, α′ = 0, filtered_vorticity = false)

Forcing due to mutual friction of a normal fluid with a Fourier-filtered superfluid velocity.

This forcing is similar to [`NormalFluidForcing`](@ref), but tries to only affect scales
within a given band `[kmin, kmax]` in Fourier space. This is achieved by a normal fluid velocity
field represented by a [`FourierBandVectorField`](@ref), and by a modified Schwarz's
equation in which only a coarse-grained superfluid flow is taken into account in the
estimation of the mutual friction term.

Concretely, the vortex line velocity according to this forcing type is:

```math
\frac{\mathrm{d}\bm{s}}{\mathrm{d}t} = \bm{v}_{\text{L}} = \bm{v}_{\text{s}} + \bm{v}_{\text{f}}
```

The forcing velocity is of the form:

```math
\bm{v}_{\text{f}} = α \bm{s}' × \bm{v}_{\text{ns}}^{>} - α' \bm{s}' × \left( \bm{s}' × \bm{v}_{\text{ns}}^{>} \right)
```

where ``\bm{v}_{\text{ns}}^{>} = \bm{v}_{\text{n}} - \bm{v}_{\text{s}}^{>}`` is a filtered slip velocity.
In practice, the filtered velocity is active within the same wavenumber range `[kmin, kmax]` where `vn` is
defined. See [`NormalFluidForcing`](@ref) for other definitions.

## Using a filtered vorticity field

To further ensure that this forcing only affects the chosen range of scales, one can pass
`filtered_vorticity = true`, which will replace the local unit tangent ``\bm{s}'`` with a
normalised coarse-grained vorticity. This corresponds to setting ``\bm{s}' = \bm{ω}^{>} / |\bm{ω}^{>}|``
where ``\bm{ω}^{>}`` is the Fourier-filtered vorticity field.
"""
struct FourierBandForcing{
        T <: AbstractFloat,
        N,  # number of dimensions (usually 3)
        VelocityField <: FourierBandVectorField{T, N},
    } <: AbstractForcing
    vn :: VelocityField
    α  :: T
    α′ :: T
    filtered_vorticity :: Bool
end

with_filtered_vorticity(f::FourierBandForcing) = f.filtered_vorticity

function FourierBandForcing(
        vn::FourierBandVectorField{T, 3}; α::Real, α′::Real = 0,
        filtered_vorticity::Bool = false,
    ) where {T <: AbstractFloat}
    FourierBandForcing(vn, T(α), T(α′), filtered_vorticity)
end

function Base.show(io::IO, f::FourierBandForcing{T}) where {T}
    (; vn, α, α′,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "FourierBandForcing{$T} with:")
    print(io, "\n$(prefix)├─ Normal velocity field: ", vn)
    print(io, "\n$(prefix)├─ Filtered superfluid vorticity: ", with_filtered_vorticity(f))
    print(io, "\n$(prefix)└─ Friction coefficients: α = ", α, " and α′ = ", α′)
end

# Here vs_grid is the superfluid velocity in Fourier space, optionally on a device (GPU).
function init_cache(f::FourierBandForcing{T, N}, cache_bs::BiotSavartCache) where {T, N}
    (; vn,) = f
    vs_grid = BiotSavart.get_longrange_field_fourier(cache_bs).field
    backend = KA.get_backend(vs_grid[1])  # CPU, CUDABackend, ROCBackend, ...
    vn_d = adapt(backend, vn)::FourierBandVectorField        # vn on the device
    vtmp_h = similar(vn)::FourierBandVectorField             # temporary buffer (on host)
    vtmp_d = adapt(backend, vtmp_h)::FourierBandVectorField  # temporary buffer (on device)
    ω_h = similar(vtmp_h)  # coarse-grained superfluid vorticity
    if !with_filtered_vorticity(f)
        empty!(ω_h)  # we don't use this field
    end
    ω_d = adapt(backend, ω_h)
    (; vn_d, vtmp_d, vtmp_h, ω_h, ω_d)
end

function update_cache!(cache, f::FourierBandForcing{T, N}, cache_bs::BiotSavartCache) where {T, N}
    (; vtmp_d, vn_d, vtmp_h, ω_h, ω_d) = cache

    vs_grid, ks_grid, σ_gaussian = let data = BiotSavart.get_longrange_field_fourier(cache_bs)
        local (; state, field, wavenumbers,) = data
        @assert state.quantity == :velocity
        field, wavenumbers, state.smoothing_scale
    end

    # (0) Copy normal fluid velocity from CPU to GPU, in case it has changed (e.g. if we're
    # using a velocity that varies in time). We don't need to do this if vn_d and f.vn are
    # the same object (typically in pure CPU simulations).
    if f.vn !== vn_d
        copyto!(vn_d, f.vn)
    end

    # (1) Copy normal fluid velocity onto buffer (device -> device copy)
    @assert vtmp_d.cs !== vn_d.cs  # coefficients are not aliased
    copyto!(vtmp_d, vn_d)

    # (2) Retrieve values from coarse-grained (Gaussian-filtered) velocity field.
    # We "unfilter" the values, similarly to the way we compute energy spectra from the long-range velocity.
    σ²_over_two = σ_gaussian^2 / 2
    @inline function op(vn, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        φ = @fastmath exp(σ²_over_two * k²)
        vs = φ * vs_filtered
        vn - vs
    end
    SyntheticFields.from_fourier_grid!(op, vtmp_d, vs_grid, ks_grid)  # vtmp_d now contains vn_d(k⃗) - vs_grid(k⃗) in [kmin, kmax]

    # Optionally compute coarse-grained vorticity.
    @inline function op_vorticity(_, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        φ = @fastmath exp(σ²_over_two * k²)
        vs = φ * vs_filtered
        im * (k⃗ × vs)
    end
    if with_filtered_vorticity(f)
        @assert length(ω_h) == length(ω_d) == length(vtmp_d)
        SyntheticFields.from_fourier_grid!(op_vorticity, ω_d, vs_grid, ks_grid)
    end

    # (3) Copy results to CPU if needed (avoided if the "device" is the CPU).
    if vtmp_h !== vtmp_d
        copyto!(vtmp_h, vtmp_d)
    end
    if ω_h !== ω_d && with_filtered_vorticity(f)
        copyto!(ω_h, ω_d)
    end

    nothing
end

function apply!(forcing::FourierBandForcing, cache, vs::AbstractVector, f::AbstractFilament; vdiff = nothing, scheduler = SerialScheduler())
    (; α, α′,) = forcing
    (; vtmp_h, ω_h,) = cache  # vtmp_h contains vₙ - vₛ at large scale
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    tforeach(eachindex(vs); scheduler) do i
        @inline
        s⃗ = f[i]
        if with_filtered_vorticity(forcing)
            ω⃗ = V(ω_h(s⃗))  # coarse-grained vorticity at s⃗
            ω_norm = sqrt(sum(abs2, ω⃗))
            # Note: in the case that the coarse-grained vorticity is exactly zero (very very
            # unlikely), we simply set s⃗′ to zero which disables the forcing.
            ω_invnorm = ifelse(iszero(ω_norm), ω_norm, 1/ω_norm)  # = 1/|ω⃗| (or zero, if |ω⃗| = 0)
            s⃗′ = (ω⃗ * ω_invnorm)::V  # replace s⃗′ with direction of coarse-grained vorticity
        else
            s⃗′ = f[i, UnitTangent()]::V
        end
        v⃗ₙₛ = V(vtmp_h(s⃗))
        if vdiff !== nothing
            vdiff[i] = v⃗ₙₛ
        end
        v_perp = s⃗′ × v⃗ₙₛ
        vf = α * v_perp
        if !iszero(α′)
            vf = vf - α′ * s⃗′ × v_perp
        end
        vs[i] = vs[i] + vf
    end
    vs
end
