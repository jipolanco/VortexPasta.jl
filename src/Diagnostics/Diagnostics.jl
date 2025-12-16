"""
    Diagnostics

Contains tools for computing different diagnostics (total energy, energy spectra, ...) from
simulation data.
"""
module Diagnostics

using ..PaddedArrays: PaddedVector
using ..Filaments:
    Filaments,
    AbstractFilament, ClosedFilament,
    Derivative, UnitTangent, CurvatureVector, Vec3,
    FilamentChunkIterator,
    knots, segments, integrate,
    number_type

using ..BiotSavart: BiotSavart, BiotSavartCache, LongRangeCache, ParamsBiotSavart, Infinity, ∞,
                    ka_generate_kernel, ka_default_workgroupsize

using Bumper: Bumper, @no_escape, @alloc
using Adapt: adapt
using KernelAbstractions:
    KernelAbstractions as KA,
    @kernel, @index, @Const, @groupsize, @localmem, @synchronize, @uniform, @ndrange
using StructArrays: StructArrays, StructArray
using LinearAlgebra: ⋅, ×

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const SingleFilamentData = AbstractVector{<:Vec3}
const SetOfFilamentsData = AbstractVector{<:SingleFilamentData}

# Trait indicating whether we can interpolate a list of values, e.g. a vector of velocities
# on filament discretisation points. A filament (or a filament-like object) can be
# interpolated, while a regular vector cannot.
struct IsInterpolable{B}
    IsInterpolable(b::Bool) = new{b}()
end

(::IsInterpolable{a} & ::IsInterpolable{b}) where {a, b} = IsInterpolable(a && b)

isinterpolable(::Type{<:AbstractVector}) = IsInterpolable(false)
isinterpolable(::Type{<:AbstractFilament}) = IsInterpolable(true)
isinterpolable(u::AbstractVector) = isinterpolable(typeof(u))  # note: this also applies to filaments (since AbstractFilament <: AbstractVector)

function _domain_volume(Ls)
    V = prod(Ls)
    if V === Infinity()
        true  # set volume to 1 for infinite domain
    else
        V
    end
end

_domain_volume(p::ParamsBiotSavart) = _domain_volume(p.Ls)

function maybe_parallelise_sum(f::F, fs::VectorOfFilaments, nthreads) where {F <: Function}
    T = number_type(fs)
    if nthreads == 1
        x = zero(T)
        for i in eachindex(fs)
            x += @inline f(i, eachindex(fs[i]))::T
        end
    else
        x_ref = Threads.Atomic{T}(zero(T))
        @sync for chunk in FilamentChunkIterator(fs; nchunks = nthreads)
            Threads.@spawn let x_local = zero(T)
                for (i, inds, _) in chunk
                    x_local += @inline f(i, inds)::T
                end
                Threads.atomic_add!(x_ref, x_local)
            end
        end
        x = x_ref[]::T
    end
    x
end

include("energy.jl")
include("energy_injection.jl")
include("helicity.jl")
include("filament_length.jl")
include("vortex_impulse.jl")
include("stretching.jl")
include("spectra.jl")
include("spectra_helicity.jl")
include("energy_flux.jl")
include("integral_scale.jl")

end
