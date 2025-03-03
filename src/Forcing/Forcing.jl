"""
    Forcing

Defines methods for injecting energy onto a system of vortices.
"""
module Forcing

using ..Filaments: AbstractFilament, UnitTangent
using ..SyntheticFields: SyntheticFields, FourierBandVectorField
using LinearAlgebra: ×

using Adapt: adapt

export AbstractForcing, NormalFluidForcing, FourierBandForcing

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
"""
function apply! end

include("normal_fluid.jl")
include("fourier_band.jl")

end
