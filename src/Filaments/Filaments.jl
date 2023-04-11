"""
    Filaments

Module for dealing with the discretisation of curves in 3D space.
"""
module Filaments

export
    ClosedFilament,
    Vec3,
    Derivative,
    nodes,
    update_coefficients!,
    normalise_derivatives,
    normalise_derivatives!,
    derivatives,
    derivative

using LinearAlgebra: norm, normalize, ⋅, ×
using StaticArrays
using StructArrays

"""
    Vec3{T}

Three-element static vector, alias for `SVector{3, T}`.

Used to describe vectors and coordinates in 3D space.
"""
const Vec3{T} = SVector{3, T}

"""
    Derivative{N}

Represents the ``N``-th order derivative operator.

Used in particular to interpolate derivatives along filaments.
"""
struct Derivative{N} end
Derivative(N::Int) = Derivative{N}()
Base.broadcastable(d::Derivative) = Ref(d)  # disable broadcasting on Derivative objects

@doc raw"""
    AbstractFilament{T} <: AbstractVector{Vec3{T}}

Abstract type representing a curve in 3D space.

The curve coordinates are parametrised as ``\bm{X}(t)`` with ``t ∈ ℝ``.

The curve is discretised by a set of *nodes* (or discretisation points)
``\bm{X}(t_i) = \bm{X}_i`` for ``i ∈ \{1, 2, …, N\}``.

See [`ClosedSplineFilament`](@ref) for a concrete implementation of `AbstractFilament`.

An `AbstractFilament` is treated as an `AbstractVector` of length `N`, in which
each element is a discretisation point `\bm{X}_i`. Therefore, one can use the
usual indexing notation to retrieve and to modify discretisation points. See
[`ClosedSplineFilament`](@ref) for some examples.
"""
abstract type AbstractFilament{T} <: AbstractVector{Vec3{T}} end

"""
    discretisation_method(f::AbstractFilament) -> DiscretisationMethod

Return the method used to discretise the filament based on its node locations.
"""
function discretisation_method end

"""
    nodes(f::AbstractFilament{T}) -> AbstractVector{T}

Return the discretisation points ``\\bm{X}_i`` of the filament.
"""
nodes(f::AbstractFilament) = f.Xs

"""
    Base.getindex(f::AbstractFilament{T}, i::Int) -> Vec3{T}

Return coordinates of discretisation point ``\\bm{X}_i``.
"""
Base.@propagate_inbounds Base.getindex(f::AbstractFilament, i::Int) = nodes(f)[i]

"""
    Base.setindex!(f::AbstractFilament{T}, v, i::Int) -> Vec3{T}

Set coordinates of discretisation point ``\\bm{X}_i``.
"""
Base.@propagate_inbounds Base.setindex!(f::AbstractFilament, v, i::Int) = nodes(f)[i] = v

Base.eltype(::Type{<:AbstractFilament{T}}) where {T} = T
Base.eltype(f::AbstractFilament) = eltype(typeof(f))
Base.size(f::AbstractFilament) = size(nodes(f))

function Base.showarg(io::IO, f::AbstractFilament, toplevel)
    toplevel || print(io, "::")
    T = eltype(f)
    disc = typeof(discretisation_method(f))
    print(io, nameof(typeof(f)), '{', T, ',', ' ', disc, '}')
end

"""
    ClosedFilament{T} <: AbstractFilament{T}

Abstract type representing a *closed* curve (a loop) in 3D space.
"""
abstract type ClosedFilament{T} <: AbstractFilament{T} end

include("discretisations.jl")
include("padded_vector.jl")

include("local/finitediff.jl")
include("local/interpolation.jl")
include("local/interp_hermite.jl")
include("local/closed_filament.jl")

include("spline/spline.jl")
include("spline/closed_filament.jl")

"""
    Filaments.init(ClosedFilament{T}, N::Integer, method::DiscretisationMethod, [args...]) -> ClosedFilament{T}

Allocate data for a closed filament with `N` discretisation points.

The element type `T` can be omitted, in which case the default `T = Float64` is used.

Depending on the type of `method`, the returned filament may be a
[`ClosedLocalFilament`](@ref) or a [`ClosedSplineFilament`](@ref).
See their respective documentations for possible optional arguments (`args...`).
"""
function init end

init(::Type{ClosedFilament}, args...) = init(ClosedFilament{Float64}, args...)

"""
    update_coefficients!(f::AbstractFilament)

Compute coefficients needed to perform inter-node interpolations and estimate
derivatives.

Uses the current locations of the filament nodes. If nodes change, this
function should be called to update the coefficients.

In the case of local Hermite interpolations, the coefficients are just the
derivatives at the discretisation points.

Note that derivatives are with respect to the (arbitrary) parametrisation
``\\bm{X}(t)``, and *not* with respect to the arclength ``ξ = ξ(t)``. In other
words, the returned derivatives do not directly correspond to the unit tangent
and curvature vectors (but they are closely related).
"""
function update_coefficients! end

"""
    normalise_derivatives(Ẋ::Vec3, Ẍ::Vec3) -> (X′, X″)
    normalise_derivatives((Ẋ, Ẍ)::NTuple)   -> (X′, X″)

Return derivatives with respect to the arc length ``ξ``, from derivatives with
respect to the parameter ``t``.

The returned derivatives satisfy:

- ``\\bm{X}' ≡ t̂`` is the **unit tangent** vector;

- ``\\bm{X}'' ≡ ρ n̂`` is the **curvature** vector, where ``n̂`` is the normal unit
  vector (with ``t̂ ⋅ n̂ = 0``) and ``ρ = R^{-1}`` is the curvature (and R the
  curvature radius).
"""
function normalise_derivatives(Ẋ::Vec3, Ẍ::Vec3)
    t̂ = normalize(Ẋ)  # unit tangent vector (= X′)
    X″ = (Ẍ - (Ẍ ⋅ t̂) * t̂) ./ sum(abs2, Ẋ)  # curvature vector
    t̂, X″
end

normalise_derivatives(derivs::NTuple{2, Vec3}) = normalise_derivatives(derivs...)

"""
    normalise_derivatives!(fil::AbstractFilament)
    normalise_derivatives!(Ẋ::AbstractVector, Ẍ::AbstractVector)

Normalise derivatives at filament nodes.

Note that filament derivatives are modified, and thus Hermite interpolations
may be incorrect after doing this. If possible, prefer using
[`normalise_derivatives`](@ref), which works on a single filament location at a time.

See [`normalise_derivatives`](@ref) for more details.
"""
function normalise_derivatives!(Ẋ::AbstractVector, Ẍ::AbstractVector)
    derivs = StructArray((Ẋ, Ẍ))
    map!(normalise_derivatives, derivs, derivs)
    (Ẋ, Ẍ)
end

normalise_derivatives!(fil::AbstractFilament) = normalise_derivatives!(derivatives(fil)...)

end
