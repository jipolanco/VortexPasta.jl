@doc raw"""
    FourierBandForcing <: AbstractForcing
    FourierBandForcing(vn::FourierBandVectorField; α)

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

where the forcing velocity is

```math
\bm{v}_{\text{f}} = α \bm{s}' × (\bm{v}_{\text{n}} - \bm{v}_{\text{s}}^{>})
```

Here ``\bm{v}_{\text{s}}^{>}`` is a Fourier (band pass) filtered superfluid velocity field.
In practice, the filtered velocity is active within the same wavenumber range `[kmin, kmax]` where `vn` is
defined.
"""
struct FourierBandForcing{
        T,
        N,  # number of dimensions (usually 3)
        VelocityField <: FourierBandVectorField{T, N},
    } <: AbstractForcing
    vn :: VelocityField
    α  :: T
end

function Base.show(io::IO, f::FourierBandForcing{T}) where {T}
    (; vn, α,) = f
    indent = get(io, :indent, 0)
    nspaces = max(indent, 1)
    spaces = " "^nspaces
    print(io, "FourierBandForcing{$T} with:")
    print(io, "\n$(spaces)├─ Normal velocity field: ", vn)
    print(io, "\n$(spaces)└─ Friction coefficient: α = ", α)
end

# Here vs_d is the superfluid velocity in Fourier space, optionally on a device (GPU).
function init_cache(f::FourierBandForcing{T, N}, vs_d::NTuple{N, AbstractArray}) where {T, N}
    (; vn,) = f
    A = eltype(vs_d)        # e.g. Array{...}, CuArray{...}, ...
    vn_d = adapt(A, vn)::FourierBandVectorField        # vn on the device
    vtmp_h = similar(vn)::FourierBandVectorField       # temporary buffer (on host)
    vtmp_d = adapt(A, vtmp_h)::FourierBandVectorField  # temporary buffer (on device)
    (; vn_d, vtmp_d, vtmp_h,)
end

function update_cache!(cache, ::FourierBandForcing{T, N}, vs_d::NTuple{N}, α_ewald::Real) where {T, N}
    (; vtmp_d, vn_d, vtmp_h,) = cache

    # (1) Copy normal fluid velocity onto buffer (device -> device copy)
    @assert vtmp_d.cs !== vn_d.cs  # coefficients are not aliased
    copyto!(vtmp_d, vn_d)

    # (2) Retrieve values from coarse-grained velocity field (Gaussian filtered with α_ewald parameter).
    # We "unfilter" the values, similarly to the way we compute energy spectra from the long-range velocity.
    @inline function op(vn, vs_filtered, k⃗)
        k² = sum(abs2, k⃗)
        β = @fastmath exp(k² / (2 * α_ewald^2))
        vs = β * vs_filtered
        vn - vs
    end
    SyntheticFields.from_fourier_grid!(op, vtmp_d, vs_d)  # vtmp_d now contains vn_d(k⃗) - vs_d(k⃗) in [kmin, kmax]

    # (3) Copy results to CPU if needed (avoided if the "device" is the CPU).
    if vtmp_h !== vtmp_d
        copyto!(vtmp_h, vtmp_d)
    end

    nothing
end

function apply!(forcing::FourierBandForcing, cache, vs::AbstractVector, f::AbstractFilament)
    (; α,) = forcing
    (; vtmp_h,) = cache  # contains vₙ - vₛ at large scale
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    for i in eachindex(vs)
        s⃗ = f[i]
        s⃗′ = f[i, UnitTangent()]
        v⃗ₙₛ = V(vtmp_h(s⃗))
        v_perp = s⃗′ × v⃗ₙₛ
        vf = α * v_perp    # velocity due to Magnus force
        vs[i] = vs[i] + vf
    end
    vs
end
