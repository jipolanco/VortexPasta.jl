"""
    Forcing

Defines methods for injecting energy onto a system of vortices.
"""
module Forcing

using ..Filaments: Filaments, AbstractFilament, UnitTangent, Derivative, UnitTangent
using ..BiotSavart: BiotSavart, BiotSavartCache
using ..SyntheticFields: SyntheticFields, FourierBandVectorField
using LinearAlgebra: ×
using KernelAbstractions: KernelAbstractions as KA
using OhMyThreads: Scheduler, SerialScheduler, tforeach, tmapreduce

using Adapt: adapt

export AbstractForcing, NormalFluidForcing, FourierBandForcing, FourierBandForcingBS

"""
    AbstractForcing

Abstract type representing a forcing method.
"""
abstract type AbstractForcing end

init_cache(forcing::Nothing, args...) = nothing  # called when forcing is disabled

"""
    Forcing.apply!(forcing::AbstractForcing, vs::AbstractVector{<:Vec3}, f::AbstractFilament; [scheduler])

Apply forcing to a single filament `f` with self-induced velocities `vs`.

At output, the `vs` vector is overwritten with the actual vortex line velocities.

The optional `scheduler` keyword can be used to parallelise computations using one of the
[schedulers defined in OhMyThreads.jl](https://juliafolds2.github.io/OhMyThreads.jl/stable/refs/api/#Schedulers).

---

    Forcing.apply!(forcing::NormalFluidForcing, vs, vn, tangents; [scheduler])

This variant can be used in the case of a [`NormalFluidForcing`](@ref) if one already has
precomputed values of the normal fluid velocity and local unit tangents at filament points.
"""
function apply! end

include("normal_fluid.jl")
include("fourier_band.jl")
include("fourier_band_bs.jl")

end
