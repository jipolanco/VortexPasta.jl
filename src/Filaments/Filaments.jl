"""
    Filaments

Module for dealing with the discretisation of curves in 3D space.
"""
module Filaments

export
    AbstractFilament,
    ClosedFilament,
    Vec3,
    Derivative,
    knots,
    knotlims,
    nodes,
    segments,
    integrate,
    update_coefficients!,
    normalise_derivatives,
    normalise_derivatives!

using ..Quadratures: AbstractQuadrature, GaussLegendreQuadrature, quadrature
using ..BasicTypes: Vec3, Derivative

using LinearAlgebra: norm, normalize, ⋅, ×
using StaticArrays
using StructArrays

# Used internally to evaluate filament coordinates or derivatives on a given
# discretisation node. This is used when calling f[i, Derivative(n)].
struct AtNode
    i :: Int
end

@doc raw"""
    AbstractFilament{T} <: AbstractVector{Vec3{T}}

Abstract type representing a curve in 3D space.

The curve coordinates are parametrised as ``\bm{X}(t)`` with ``t ∈ ℝ``.

The curve is discretised by a set of *nodes* (or discretisation points)
``\bm{X}(t_i) = \bm{X}_i`` for ``i ∈ \{1, 2, …, N\}``.

See [`ClosedSplineFilament`](@ref) for a concrete implementation of `AbstractFilament`.

An `AbstractFilament` is treated as an `AbstractVector` of length `N`, in which
each element is a discretisation point ``\bm{X}_i``. Therefore, one can use the
usual indexing notation to retrieve and to modify discretisation points. See
[`ClosedSplineFilament`](@ref) for some examples.

# Extended help

## Evaluating coordinates and derivatives

Different ways are proposed of evaluating filament coordinates and derivatives
along a filament `f`, depending on whether one wants values on discretisation
points or in-between them.

In short, square brackets `f[...]` should be used to evaluate on filament nodes,
while round brackets `f(...)` to evaluate in-between nodes.

### Values on discretisation points

Coordinates ``\bm{X}_i`` of discretisation points can be simply obtained by
indexing the filament object:

    X = f[i]

Derivatives at discretisation points can be similarly obtained by doing:

    X′ = f[i, Derivative(1)]
    X″ = f[i, Derivative(2)]

(Note that this also works with `Derivative(0)`, in which case it's the same as `f[i]`.)

### Values in-between discretisation points

In this case, one wants to evaluate a value in-between two discretisation
points, i.e. for ``t ∈ [t_i, t_{i + 1}]`` for some index `i`.
Two options are proposed:

- if one knows the parametrisation of the filament (see
  [`knots`](@ref) to obtain the parametrisation knots), then one can evaluate
  using a value of `t`:

      X  = f(t)
      X′ = f(t, Derivative(1))
      X″ = f(t, Derivative(2))

- alternatively, if one knows the index `i` associated to the segment of
  interest, then one can do

      X  = f(i, ζ)
      X′ = f(i, ζ, Derivative(1))
      X″ = f(i, ζ, Derivative(2))

  where `ζ` must be in ``[0, 1]``, and the two limits correspond to knots ``t_i`` and ``t_{i + 1}``.
  This is convenient if one wants to evaluate, say, right in the middle between
  two discretisation points, in which case one would choose `ζ = 0.5`.

!!! note "Derivatives"

    In all cases, derivatives are computed with respect to the parametrisation ``t``.
    In particular, in `f(i, ζ, Derivative(1))`, the derivative is *not* with respect to ``ζ``.

!!! note "Derivative normalisation"

    Since ``t`` is a rough approximation for the arc length ``ξ``, first
    derivatives almost represent the **unit tangent vector** to the filament, and
    second derivatives are a rough approximation of the local **curvature vector**.

    One should use [`normalise_derivatives`](@ref) if one wants a
    more accurate estimation, which takes into account the differences between
    ``t`` and ``ξ``.
"""
abstract type AbstractFilament{T} <: AbstractVector{Vec3{T}} end

# This is needed since eltype(f) == Vec3{T}
Base.similar(f::AbstractFilament, ::Type{Vec3{T}}, dims::Dims{1}) where {T} =
    similar(f, T, dims)

"""
    ClosedFilament{T} <: AbstractFilament{T}

Abstract type representing a *closed* curve (a loop) in 3D space.
"""
abstract type ClosedFilament{T} <: AbstractFilament{T} end

"""
    discretisation_method(f::AbstractFilament) -> DiscretisationMethod

Return the method used to discretise the filament based on its node locations.
"""
function discretisation_method end

"""
    nodes(f::AbstractFilament{T}) -> AbstractVector{T}

Return the nodes (or discretisation points) ``\\bm{X}_i`` of the filament.
"""
nodes(f::AbstractFilament) = f.Xs

"""
    knots(f::AbstractFilament{T}) -> AbstractVector{T}

Return parametrisation knots ``t_i`` of the filament.
"""
knots(f::AbstractFilament) = f.ts

"""
    knotlims(f::AbstractFilament) -> (t_begin, t_end)

Return limits within which the filament can be evaluated.
"""
function knotlims end

function knotlims(f::ClosedFilament)
    ts = knots(f) :: PaddedVector
    ts[begin], ts[end + 1]  # we can index at end+1 thanks to padding
end

"""
    Base.getindex(f::AbstractFilament{T}, i::Int, [Derivative(n)]) -> Vec3{T}

Return coordinates of discretisation point ``\\bm{X}_i``.

One may also obtain derivatives at point ``\\bm{X}_i`` by passing an optional
[`Derivative`](@ref).
"""
Base.@propagate_inbounds Base.getindex(f::AbstractFilament, i::Int) = nodes(f)[i]

Base.@propagate_inbounds Base.getindex(f::AbstractFilament, i::Int, d::Derivative) =
    f(AtNode(i), d)

"""
    Base.setindex!(f::AbstractFilament{T}, v, i::Int) -> Vec3{T}

Set coordinates of discretisation point ``\\bm{X}_i``.
"""
Base.@propagate_inbounds Base.setindex!(f::AbstractFilament, v, i::Int) = nodes(f)[i] = v

Base.eltype(::Type{<:AbstractFilament{T}}) where {T} = Vec3{T}  # type returned when indexing into a filament
Base.eltype(f::AbstractFilament) = eltype(typeof(f))
Base.size(f::AbstractFilament) = size(nodes(f))

function Base.showarg(io::IO, f::AbstractFilament, toplevel)
    toplevel || print(io, "::")
    T = eltype(f)
    disc = typeof(discretisation_method(f))
    print(io, nameof(typeof(f)), '{', T, ',', ' ', disc, '}')
end

include("discretisations.jl")
include("padded_vector.jl")
include("segments.jl")
include("integrate.jl")

include("local/finitediff.jl")
include("local/interpolation.jl")
include("local/interp_hermite.jl")
include("local/closed_filament.jl")

include("spline/spline.jl")
include("spline/closed_filament.jl")

include("makie_recipes.jl")

"""
    Filaments.init(ClosedFilament{T}, N::Integer, method::DiscretisationMethod, [args...]) -> ClosedFilament{T}

Allocate data for a closed filament with `N` discretisation points.

The element type `T` can be omitted, in which case the default `T = Float64` is used.

Depending on the type of `method`, the returned filament may be a
[`ClosedLocalFilament`](@ref) or a [`ClosedSplineFilament`](@ref).
See their respective documentations for possible optional arguments (`args...`).
"""
init(::Type{ClosedFilament}, N::Integer, args...) = init(ClosedFilament{Float64}, N, args...)

"""
    Filaments.init(ClosedFilament, points::AbstractVector{<:Vec3}, method::DiscretisationMethod, [args...])

Initialise new filament with the chosen discretisation points.

Note that [`update_coefficients!`](@ref) does not need to be called after using
this variant (until node locations change, of course!).
"""
function init(::Type{ClosedFilament}, positions::AbstractVector{<:Vec3}, args...)
    T = eltype(eltype(positions))
    f = init(ClosedFilament{T}, length(positions), args...)
    copy!(nodes(f), positions)
    update_coefficients!(f)
    f
end

"""
    update_coefficients!(f::AbstractFilament)

Compute coefficients needed to perform inter-node interpolations and estimate
derivatives.

Uses the current locations of the filament nodes. If nodes change, this
function should be called to update the coefficients.

In the case of local Hermite interpolations, the coefficients are just the
derivatives at the discretisation points.

Note that derivatives are with respect to the (arbitrary) parametrisation
``\\bm{X}(t)``, and *not* with respect to the arc length ``ξ = ξ(t)``. In other
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
    normalise_derivatives!(Ẋ::AbstractVector, Ẍ::AbstractVector)

Normalise vectors containing derivatives at filament locations.

If possible, prefer using [`normalise_derivatives`](@ref), which works on a
single filament location at a time.

See [`normalise_derivatives`](@ref) for more details.
"""
function normalise_derivatives!(Ẋ::AbstractVector, Ẍ::AbstractVector)
    derivs = StructArray((Ẋ, Ẍ))
    map!(normalise_derivatives, derivs, derivs)
    (Ẋ, Ẍ)
end

# Update filament parametrisation knots `ts` from node coordinates `Xs`.
function _update_knots_periodic!(ts::PaddedVector, Xs::PaddedVector)
    @assert eachindex(ts) == eachindex(Xs)
    ts[begin] = 0
    inds = eachindex(ts)
    @assert npad(ts) == npad(Xs) ≥ 1
    @inbounds for i ∈ inds
        ts[i + 1] = ts[i] + norm(Xs[i + 1] - Xs[i])
    end
    L = ts[end + 1] - ts[begin]  # knot period
    pad_periodic!(ts, L)
end

# TESTING / EXPERIMENTAL
# This function may be removed in the future.
# The idea is to update the parametrisation of `f` to follow more closely the
# actual arc lengths of the filament. Not sure if it's worth it...
function recompute_parametrisation!(f::ClosedFilament)
    m = interpolation_method(f)
    _recompute_parametrisation!(m, f)
end

# In the case of straight segments (linear interpolation), the parametrisation
# cannot be improved from its initial estimation.
_recompute_parametrisation!(::HermiteInterpolation{0}, f::AbstractFilament) = f

function _recompute_parametrisation!(::Any, f::AbstractFilament)
    (; ts,) = f
    quad = GaussLegendreQuadrature(4)
    @assert npad(ts) ≥ 1
    tnext = ts[begin]
    for i ∈ eachindex(ts)
        # Estimate arc length from ts[i] to ts[i + 1]
        ℓ = integrate(f, i, quad) do ζ
            norm(f(i, ζ, Derivative(1)))  # = ∂X/∂t
        end
        @assert ℓ ≥ ts[i + 1] - ts[i]
        ts[i] = tnext
        tnext += ℓ  # this will be the new value of ts[i + 1], but we can't update it yet...
    end
    L = tnext - ts[begin]  # full length of the filament
    pad_periodic!(ts, L)
    _update_coefficients_only!(f)
    f
end

end
