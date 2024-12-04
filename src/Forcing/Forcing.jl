"""
    Forcing

Defines methods for injecting energy onto a system of vortices.
"""
module Forcing

using ..Filaments: AbstractFilament, UnitTangent
using ..SyntheticFields: SyntheticFields, SyntheticVectorField
using LinearAlgebra: ×

export AbstractForcing, NormalFluidForcing

"""
    AbstractForcing

Abstract type representing a forcing method.
"""
abstract type AbstractForcing end

"""
    Forcing.apply!(forcing::AbstractForcing, vs::AbstractVector, f::AbstractFilament)

Apply forcing to a single filament `f`, modifying vortex velocities `vs`.
"""
function apply! end

@doc raw"""
    NormalFluidForcing <: AbstractForcing
    NormalFluidForcing(field::SyntheticVectorField; α, α′ = 0)

Forcing due to mutual friction with a normal fluid.

The normal fluid is represented by a synthetic velocity field ``\bm{v}_{\text{n}}(\bm{x}, \bm{t})``
from the [`SyntheticFields`](@ref) module.

This type of forcing defines an external velocity ``\bm{v}_{\text{f}}`` affecting vortex
motion, so that the effective vortex velocity is

```math
\frac{\mathrm{d}\bm{s}}{\mathrm{d}t} = \bm{v}_{\text{s}} + \bm{v}_{\text{f}}.
```

Here ``\bm{v}_{\text{s}}`` is the self-induced velocity obtained by applying Biot–Savart's law.

The forcing velocity is of the form:

```math
\bm{v}_{\text{f}} =
α \bm{s}' × (\bm{v}_{\text{n}} - \bm{v}_{\text{s}})
- α' \bm{s}' × \left[ \bm{s}' × (\bm{v}_{\text{n}} - \bm{v}_{\text{s}}) \right]
```

where ``\bm{s}'`` is the local unit tangent vector, and ``α`` and ``α'`` are non-dimensional
coefficients representing the intensity of Magnus and drag forces.

# Example

Define a mutual friction forcing based on a large-scale normal fluid velocity field
(see [`SyntheticFields.FourierBandVectorField`](@ref)):

```jldoctest
julia> using VortexPasta.Forcing: NormalFluidForcing

julia> using VortexPasta.SyntheticFields: SyntheticFields, FourierBandVectorField

julia> using Random: Xoshiro

julia> rng = Xoshiro(42);  # initialise random number generator (optional, but recommended)

julia> Ls = (2π, 2π, 2π);  # domain dimensions

julia> vn_rms = 1.0;  # typical magnitude (rms value) of normal fluid velocity components

julia> vn = FourierBandVectorField(undef, Ls; kmin = 0.1, kmax = 1.5)  # create field with non-zero Fourier wavevectors kmin ≤ |k⃗| ≤ kmax
FourierBandVectorField{Float64, 3} with 9 independent Fourier coefficients in |k⃗| ∈ [1.0, 1.4142]

julia> SyntheticFields.init_coefficients!(rng, vn, vn_rms);  # randomly set non-zero Fourier coefficients of the velocity field

julia> forcing = NormalFluidForcing(vn; α = 0.8, α′ = 0)
NormalFluidForcing{Float64, 3} with:
 - Magnus force coefficient: α = 0.8
 - Drag force coefficient: α′ = 0.0
 - Normal velocity field: FourierBandVectorField{Float64, 3} with 9 independent Fourier coefficients in |k⃗| ∈ [1.0, 1.4142]
```

"""
struct NormalFluidForcing{
        T <: Real, N,
        VelocityField <: SyntheticVectorField{T, N}
    } <: AbstractForcing
    vn :: VelocityField
    α  :: T
    α′ :: T
end

function NormalFluidForcing(vn::SyntheticVectorField{T}; α::Real, α′::Real = 0) where {T}
    NormalFluidForcing(vn, T(α), T(α′))
end

function Base.show(io::IO, f::NormalFluidForcing{T, N}) where {T, N}
    (; vn, α, α′,) = f
    indent = get(io, :indent, 0)
    nspaces = max(indent, 1)
    spaces = " "^nspaces
    print(io, "NormalFluidForcing{$T, $N} with:")
    print(io, "\n$(spaces)├─ Magnus force coefficient: α = ", α)
    print(io, "\n$(spaces)├─ Drag force coefficient: α′ = ", α′)
    print(io, "\n$(spaces)└─ Normal velocity field: ", vn)
end

function apply!(forcing::NormalFluidForcing, vs::AbstractVector, f::AbstractFilament)
    length(vs) == length(f) || throw(DimensionMismatch("length of filament and velocity vector don't match"))
    (; vn, α, α′,) = forcing
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    # TODO: parallelise?
    for i ∈ eachindex(vs, f)
        vs_i = vs[i]  # vortex velocity before mutual friction forcing
        s⃗ = f[i]                  # vortex point location
        s⃗′ = f[i, UnitTangent()]  # local tangent
        vn_i = V(vn(s⃗))        # evaluate normal fluid velocity at s⃗
        v_ns = vn_i - vs_i     # slip velocity
        v_perp = s⃗′ × v_ns
        vf = α * v_perp      # velocity due to Magnus force
        if !iszero(α′)
            vf = vf - α′ * s⃗′ × v_perp  # velocity due to drag force (it's quite common to set α′ = 0)
        end
        vs[i] = vs[i] + vf
    end
    vs
end

end
