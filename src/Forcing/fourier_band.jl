@doc raw"""
    FourierBandForcing <: AbstractForcing
    FourierBandForcing(vn::FourierBandVectorField; α, α′ = 0, filtered_vorticity = false)

Forcing due to _large scale_ mutual friction with a normal fluid.

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
\bm{v}_{\text{f}} = α \bm{s}' × \bm{v}_{\text{ns}}^{>} - α′ \bm{s}' × \left[ \bm{s}' × \bm{v}_{\text{ns}}^{>} \right]
```

where ``\bm{v}_{\text{ns}}^{>} = \bm{v}_{\text{n}} - \bm{v}_{\text{s}}^{>}`` is a coarse-grained slip velocity.
In practice, the coarse-grained velocity is active within the same wavenumber range `[kmin, kmax]` where `vn` is
defined. See [`NormalFluidForcing`](@ref) for other definitions.

## Using a filtered vorticity field

To further ensure that this forcing only affects the large scales
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
    indent = get(io, :indent, 0)
    nspaces = max(indent, 1)
    spaces = " "^nspaces
    print(io, "FourierBandForcing{$T} with:")
    print(io, "\n$(spaces)├─ Normal velocity field: ", vn)
    print(io, "\n$(spaces)├─ Filtered superfluid vorticity: ", with_filtered_vorticity(f))
    print(io, "\n$(spaces)└─ Friction coefficients: α = ", α, " and α′ = ", α′)
end

# Here vs_d is the superfluid velocity in Fourier space, optionally on a device (GPU).
function init_cache(f::FourierBandForcing{T, N}, cache_bs::BiotSavartCache) where {T, N}
    (; vn,) = f
    vs_d = BiotSavart.get_longrange_field_fourier(cache_bs).field
    A = eltype(vs_d).name.wrapper  # e.g. Array, CuArray, ... (XXX: relies on Julia internals!)
    vn_d = adapt(A, vn)::FourierBandVectorField        # vn on the device
    vtmp_h = similar(vn)::FourierBandVectorField       # temporary buffer (on host)
    vtmp_d = adapt(A, vtmp_h)::FourierBandVectorField  # temporary buffer (on device)
    ω_h = similar(vtmp_h)  # coarse-grained superfluid vorticity
    if !with_filtered_vorticity(f)
        empty!(ω_h)  # we don't use this field
    end
    ω_d = adapt(A, ω_h)
    (; vn_d, vtmp_d, vtmp_h, ω_h, ω_d)
end

function update_cache!(cache, f::FourierBandForcing{T, N}, cache_bs::BiotSavartCache) where {T, N}
    (; vtmp_d, vn_d, vtmp_h, ω_h, ω_d) = cache

    vs_d, ks_grid = let data = BiotSavart.get_longrange_field_fourier(cache_bs)
        local (; state, field, wavenumbers,) = data
        @assert state.quantity == :velocity
        @assert state.smoothed == true
        field, wavenumbers
    end
    α_ewald = cache_bs.params.α

    # (1) Copy normal fluid velocity onto buffer (device -> device copy)
    @assert vtmp_d.cs !== vn_d.cs  # coefficients are not aliased
    copyto!(vtmp_d, vn_d)

    # (2) Retrieve values from coarse-grained velocity field (Gaussian filtered with α_ewald parameter).
    # We "unfilter" the values, similarly to the way we compute energy spectra from the long-range velocity.
    inv_four_α² = 1 / (4 * α_ewald * α_ewald)
    @inline function op(vn, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        φ = @fastmath exp(k² * inv_four_α²)
        vs = φ * vs_filtered
        vn - vs
    end
    SyntheticFields.from_fourier_grid!(op, vtmp_d, vs_d, ks_grid)  # vtmp_d now contains vn_d(k⃗) - vs_d(k⃗) in [kmin, kmax]

    # Optionally compute coarse-grained vorticity.
    @inline function op_vorticity(_, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        φ = @fastmath exp(k² * inv_four_α²)
        vs = φ * vs_filtered
        im * (k⃗ × vs)
    end
    if with_filtered_vorticity(f)
        @assert length(ω_h) == length(ω_d) == length(vtmp_d)
        SyntheticFields.from_fourier_grid!(op_vorticity, ω_d, vs_d, ks_grid)
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

function apply!(forcing::FourierBandForcing, cache, vs::AbstractVector, f::AbstractFilament)
    (; α, α′,) = forcing
    (; vtmp_h, ω_h,) = cache  # contains vₙ - vₛ at large scale
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    for i in eachindex(vs)
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
        v_perp = s⃗′ × v⃗ₙₛ
        vf = α * v_perp
        if !iszero(α′)
            vf = vf - α′ * s⃗′ × v_perp
        end
        vs[i] = vs[i] + vf
    end
    vs
end
