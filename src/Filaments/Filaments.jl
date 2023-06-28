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

    # Geometric quantities
    UnitTangent,
    CurvatureVector,
    CurvatureScalar,
    CurvatureBinormal,

    # Reconnection criteria
    NoReconnections,
    ReconnectBasedOnDistance,

    knots,
    knotlims,
    end_to_end_offset,
    minimum_knot_increment,
    maximum_knot_increment,
    nodes,
    segments,
    integrate,
    update_coefficients!,
    normalise_derivatives,
    normalise_derivatives!

using ..Quadratures: AbstractQuadrature, GaussLegendre, quadrature
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

For convenience, other geometric quantities can be evaluated in a similar way:

    ρ⃗ = f[i, CurvatureVector()]

(See [`GeometricQuantity`](@ref) for available quantities.)

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

For convenience, other geometric quantities can be evaluated in a similar way:

    ρ⃗ = f(t, CurvatureVector())
    ρ⃗ = f(i, ζ, CurvatureVector())

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

function Base.copyto!(v::F, u::F) where {F <: AbstractFilament}
    map(copyto!, allvectors(v), allvectors(u))
    v
end

function Base.resize!(f::AbstractFilament, n::Integer)
    map(v -> resize!(v, n), allvectors(f))
    f
end

function Base.sizehint!(f::AbstractFilament, n::Integer)
    map(v -> sizehint!(v, n), allvectors(f))
    f
end

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
    minimum_knot_increment(f::AbstractFilament) -> Real
    minimum_knot_increment(fs::AbstractVector{<:AbstractFilament}) -> Real

Return the minimum increment ``Δt = t_{i + 1} - t_{i}`` between filament knots.

The second form allows to estimate the minimum increment among a vector of filaments.

This is generally a good approximation for the minimum segment length.
"""
minimum_knot_increment(f::AbstractFilament) = reduce_knot_increments(min, f)
minimum_knot_increment(fs::AbstractVector{<:AbstractFilament}) = minimum(minimum_knot_increment, fs)

"""
    maximum_knot_increment(f::AbstractFilament) -> Real
    maximum_knot_increment(fs::AbstractVector{<:AbstractFilament}) -> Real

Return the maximum increment ``Δt = t_{i + 1} - t_{i}`` between filament knots.

The second form allows to estimate the maximum increment among a vector of filaments.

This is generally a good approximation for the maximum segment length.
"""
maximum_knot_increment(f::AbstractFilament) = reduce_knot_increments(max, f)
maximum_knot_increment(fs::AbstractVector{<:AbstractFilament}) = maximum(maximum_knot_increment, fs)

function reduce_knot_increments(g::G, f::AbstractFilament) where {G <: Function}
    ts = knots(f)
    mapreduce(g, eachindex(segments(f))) do i
        @inbounds ts[i + 1] - ts[i]
    end
end

"""
    knotlims(f::AbstractFilament) -> (t_begin, t_end)

Return limits within which the filament can be evaluated.
"""
function knotlims end

function knotlims(f::ClosedFilament)
    ts = knots(f) :: PaddedVector
    ts[begin], ts[end + 1]  # we can index at end+1 thanks to padding
end

Base.checkbounds(::Type{Bool}, f::AbstractFilament, I...) = checkbounds(Bool, nodes(f), I...)

"""
    end_to_end_offset(f::ClosedFilament{T}) -> Vec3{T}

Return the end-to-end offset `Δ⃗ = f[end + 1] - f[begin]` of a "closed" filament.

For actually closed filaments, the end-to-end offset is zero. However, `ClosedFilament` also
supports the case of infinite (but unclosed) filaments, which infinitely extend along one or
more Cartesian directions. The restriction imposed by `ClosedFilament` is that infinite
filaments repeat themselves periodically, such that `f[i + m * N] == f[i] + m * Δ⃗` where `N`
is the `length` of the filament (i.e. the number of degrees of freedom, or the total number
of *independent* filament nodes).
"""
end_to_end_offset(f::ClosedFilament) = f.Xoffset

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
include("padded_array.jl")
include("segments.jl")
include("integrate.jl")
include("min_distance.jl")
include("refinement.jl")
include("reconnections.jl")

include("local/interpolation.jl")
include("local/interp_hermite.jl")
include("local/finitediff.jl")
include("local/closed_filament.jl")

include("spline/spline.jl")
include("spline/closed_filament.jl")

include("quantities.jl")
include("utils.jl")

include("makie_recipes.jl")

"""
    Base.getindex(f::AbstractFilament{T}, i::Int) -> Vec3{T}
    Base.getindex(f::AbstractFilament{T}, i::Int, ::Derivative{n}) -> Vec3{T}
    Base.getindex(f::AbstractFilament{T}, i::Int, ::GeometricQuantity)

Return coordinates of discretisation point ``\\bm{X}_i``.

One may also obtain derivatives and other geometric quantities at point ``\\bm{X}_i``
by passing an optional [`Derivative`](@ref) or [`GeometricQuantity`](@ref).
"""
Base.@propagate_inbounds Base.getindex(f::AbstractFilament, i::Int) = nodes(f)[i]

Base.@propagate_inbounds Base.getindex(
    f::AbstractFilament, i::Int, d::Union{Derivative, GeometricQuantity},
) = f(AtNode(i), d)

"""
    Filaments.init(
        ClosedFilament{T}, N::Integer, method::DiscretisationMethod, [args...];
        offset = zero(Vec3{T}),
    ) -> ClosedFilament{T}

Allocate data for a closed filament with `N` discretisation points.

The element type `T` can be omitted, in which case the default `T = Float64` is used.

The optional `offset` keyword argument allows to define a filament with a spatial
offset between points `f[i]` and `f[i + N]`. By default the offset is zero,
meaning that the filament is a closed loop. This can be used for defining
infinite (so not really closed) filaments living in periodic domains.

Depending on the type of `method`, the returned filament may be a
[`ClosedLocalFilament`](@ref) or a [`ClosedSplineFilament`](@ref).
See their respective documentations for possible optional arguments (`args...`).
"""
init(::Type{ClosedFilament}, N::Integer, args...; kws...) =
    init(ClosedFilament{Float64}, N, args...; kws...)

"""
    Filaments.init(
        ClosedFilament, points::AbstractVector{<:Vec3}, method::DiscretisationMethod,
        [args...];
        [kws...],
    )

Initialise new filament with the chosen discretisation points.

Note that [`update_coefficients!`](@ref) does not need to be called after using
this variant (until node locations change, of course!).
"""
function init(::Type{ClosedFilament}, positions::AbstractVector{<:Vec3}, args...; kws...)
    T = eltype(eltype(positions))
    f = init(ClosedFilament{T}, length(positions), args...; kws...)
    copy!(nodes(f), positions)
    update_coefficients!(f)
    f
end

"""
    change_offset(f::ClosedFilament, offset::Vec3{<:Real}) -> f′

Change spatial offset between a filament start and endpoints.

See [`init`](@ref) for more details on optional spatial offsets.

This function is allocation-free. It returns a new filament which shares the same
arrays as `f`, and only differs in the offset. Modifying nodes of the returned
filament also modifies nodes of `f`.
"""
function change_offset end

"""
    update_coefficients!(f::AbstractFilament; knots = nothing)

Compute coefficients needed to perform inter-node interpolations and estimate
derivatives.

Uses the current locations of the filament nodes. If nodes change, this
function should be called to update the coefficients.

By default, this function also updates the parametrisation knots ``t_i``
according to the current node positions. One can override this by passing a
`knots` vector as a keyword argument.

This function will fail if the number of filament nodes is smaller than that
required by the discretisation method. For instance, closed filaments
discretised using cubic splines must have at least 3 nodes.
"""
function update_coefficients! end

function update_coefficients!(f::ClosedFilament; knots = nothing)
    (; ts, Xs,) = f
    Xoffset = end_to_end_offset(f)
    M = npad(Xs)
    check_nodes(f)

    # 1. Periodically pad Xs.
    pad_periodic!(Xs, Xoffset)

    # 2. Compute parametrisation knots `ts`.
    @assert M == npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts
    _update_knots_periodic!(ts, Xs, knots)

    # 3. Estimate coefficients needed for derivatives and interpolations.
    _update_coefficients_only!(f)

    f
end

"""
    check_nodes(Bool, f::AbstractFilament) -> Bool
    check_nodes(f::AbstractFilament)

Check whether current filament nodes are compatible with the filament discretisation method.

In its first form, this function returns `false` in case of incompatibility,
while it throws an error in its second form.

For now, the only requirement is that the number of nodes must be larger than
some small value. In particular, one can't have a closed filament with less
than 3 nodes (but the specific discretisation method might impose some other small value).
"""
check_nodes(::Type{Bool}, f::AbstractFilament) = _check_nodes(Bool, nodes(f))
check_nodes(f::AbstractFilament) = _check_nodes(Nothing, nodes(f))

function _check_nodes(::Type{T}, Xs::PaddedVector) where {T}
    M = max(npad(Xs), 3)  # can't have a closed filament with less than 3 nodes!
    if length(Xs) < M
        if T === Bool
            return false
        else
            error(lazy"number of nodes in filament ($(length(Xs))) is below the allowed minimum ($M)")
        end
    end
    true
end

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
function _update_knots_periodic!(ts::PaddedVector, Xs::PaddedVector, ::Nothing = nothing)
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

# In this case, override the computation of knots and copy the input knot vector.
_update_knots_periodic!(ts::PaddedVector, Xs::PaddedVector, ts_in::AbstractVector) =
    copyto!(ts, ts_in)

end
