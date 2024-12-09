"""
    Forcing

Defines methods for injecting energy onto a system of vortices.
"""
module Forcing

using ..Filaments: AbstractFilament, UnitTangent
using ..SyntheticFields: SyntheticFields  # for docs only
using LinearAlgebra: ×

export AbstractForcing, NormalFluidForcing

"""
    AbstractForcing

Abstract type representing a forcing method.
"""
abstract type AbstractForcing end

"""
    Forcing.apply!(forcing::AbstractForcing, vs::AbstractVector{<:Vec3}, f::AbstractFilament)

Apply forcing to a single filament `f` with self-induced velocities `vs`.

At output, the `vs` vector is overwritten with the actual vortex line velocities.

---

    Forcing.apply!(forcing::NormalFluidForcing, vs, vn, tangents)

This variant can be used in the case of a [`NormalFluidForcing`](@ref) if one already has
precomputed values of the normal fluid velocity and local unit tangents at filament points.

The normal fluid velocities may be computed using [`Forcing.get_velocities!`](@ref).
"""
function apply! end

@doc raw"""
    NormalFluidForcing <: AbstractForcing
    NormalFluidForcing(vn::Function; α, α′ = 0)

Forcing due to mutual friction with a normal fluid.

The normal fluid is represented by a function `vn` which should take a position
`x⃗::SVector{N, T}` and return a velocity `v⃗::SVector{N, T}` (`N` is the number of
dimensions, usually `N = 3`).

In particular, the function could be a synthetic velocity field from the
[`SyntheticFields`](@ref) module (see below for examples).

This type of forcing defines an external velocity ``\bm{v}_{\text{f}}`` affecting vortex
motion, so that the actual vortex velocity ``\bm{v}_{\text{L}}`` is

```math
\frac{\mathrm{d}\bm{s}}{\mathrm{d}t} = \bm{v}_{\text{L}} = \bm{v}_{\text{s}} + \bm{v}_{\text{f}}.
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
NormalFluidForcing{Float64} with:
 ├─ Magnus force coefficient: α = 0.8
 ├─ Drag force coefficient: α′ = 0.0
 └─ Normal velocity field: FourierBandVectorField{Float64, 3} with 9 independent Fourier coefficients in |k⃗| ∈ [1.0, 1.4142]
```

"""
struct NormalFluidForcing{
        T,
        VelocityField <: Function,  # should return an SVector{N, T}
    } <: AbstractForcing
    vn :: VelocityField
    α  :: T
    α′ :: T
end

function NormalFluidForcing(vn::F; α::T, α′::Real = 0) where {T <: AbstractFloat, F <: Function}
    NormalFluidForcing(vn, T(α), T(α′))
end

function Base.show(io::IO, f::NormalFluidForcing{T}) where {T}
    (; vn, α, α′,) = f
    indent = get(io, :indent, 0)
    nspaces = max(indent, 1)
    spaces = " "^nspaces
    print(io, "NormalFluidForcing{$T} with:")
    print(io, "\n$(spaces)├─ Magnus force coefficient: α = ", α)
    print(io, "\n$(spaces)├─ Drag force coefficient: α′ = ", α′)
    print(io, "\n$(spaces)└─ Normal velocity field: ", vn)
end

"""
    Forcing.get_velocities!(forcing::NormalFluidForcing, vn::AbstractVector{<:Vec3}, f::AbstractFilament)

Evaluate normal fluid velocity at filament nodes.

Results are written to `vn`, which should have the same length as the filament `f`.
"""
function get_velocities!(forcing::NormalFluidForcing, vn::AbstractVector, f::AbstractFilament)
    for i ∈ eachindex(vn, f)
        vn[i] = forcing.vn(f[i])
    end
    vn
end

function apply!(forcing::NormalFluidForcing, vs::AbstractVector, f::AbstractFilament)
    eachindex(vs) == eachindex(f) || throw(DimensionMismatch("lengths of filament and velocity vectors don't match"))
    @inline get_at_node(i) = (v⃗ₙ = forcing.vn(f[i]), s⃗′ = f[i, UnitTangent()])
    _apply!(get_at_node, forcing, vs)
end

function apply!(forcing::NormalFluidForcing, vs::AbstractVector, vn::AbstractVector, tangents::AbstractVector)
    eachindex(vs) == eachindex(vn) == eachindex(tangents) || throw(DimensionMismatch("lengths of vectors don't match"))
    @inline get_at_node(i) = @inbounds (v⃗ₙ = vn[i], s⃗′ = tangents[i],)
    _apply!(get_at_node, forcing, vs)
end

# The first argument is a `get_at_node(i)` function which returns a NamedTuple with fields
# (v⃗ₛ, v⃗ₙ, s⃗′) with the superfluid velocity, normal fluid velocity and the local unit
# tangent at the node `i` of a filament. This is useful if those quantities have been
# precomputed.
function _apply!(get_at_node::F, forcing::NormalFluidForcing, vs::AbstractVector) where {F <: Function}
    (; α, α′,) = forcing
    V = eltype(vs)  # usually Vec3{T} = SVector{3, T}
    for i ∈ eachindex(vs)
        (; v⃗ₙ, s⃗′,) = @inline get_at_node(i)  # superfluid velocity, normal fluid velocity and unit tangent
        v⃗ₛ = vs[i]
        vₙₛ = V(v⃗ₙ) - V(v⃗ₛ)   # slip velocity
        v_perp = s⃗′ × vₙₛ
        vf = α * v_perp    # velocity due to Magnus force
        if !iszero(α′)
            vf = vf - α′ * s⃗′ × v_perp  # velocity due to drag force (it's quite common to set α′ = 0)
        end
        vs[i] = v⃗ₛ + vf
    end
    vs
end

end
