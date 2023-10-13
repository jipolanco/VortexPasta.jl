"""
    Filaments

Module for dealing with the discretisation of curves in 3D space.
"""
module Filaments

export
    AbstractFilament,
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
    redistribute_nodes!,
    normalise_derivatives,
    normalise_derivatives!

using ..Quadratures: Quadratures, AbstractQuadrature, GaussLegendre, NoQuadrature, quadrature
using ..BasicTypes: Vec3, Derivative
using ..PredefinedCurves: PredefinedCurves  # just for documentation

# Load PaddedVector and some non-exported functions associated to PaddedVector
using ..PaddedArrays: PaddedVector, pad_periodic!, FromCentre, FromRight
import ..PaddedArrays: npad  # overloaded for DiscretisationMethod

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

The curve coordinates are parametrised as ``\bm{s}(t)`` with ``t ∈ ℝ``.

The curve is discretised by a set of *nodes* (or discretisation points)
``\bm{s}(t_i) = \bm{s}_i`` for ``i ∈ \{1, 2, …, N\}``.

See [`ClosedFilament`](@ref) for a concrete implementation of `AbstractFilament`.

An `AbstractFilament` is treated as an `AbstractVector` of length `N`, in which
each element is a discretisation point ``\bm{s}_i``. Therefore, one can use the
usual indexing notation to retrieve and to modify discretisation points. See
[`ClosedFilament`](@ref) for some examples.

# Extended help

## Evaluating coordinates and derivatives

Different ways are proposed of evaluating filament coordinates and derivatives
along a filament `f`, depending on whether one wants values on discretisation
points or in-between them.

In short, square brackets `f[...]` should be used to evaluate on filament nodes,
while round brackets `f(...)` to evaluate in-between nodes.

### Values on discretisation points

Coordinates ``\bm{s}_i`` of discretisation points can be simply obtained by
indexing the filament object:

    s⃗ = f[i]

Derivatives at discretisation points can be similarly obtained by doing:

    s⃗′ = f[i, Derivative(1)]
    s⃗″ = f[i, Derivative(2)]

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

      s⃗  = f(t)
      s⃗′ = f(t, Derivative(1))
      s⃗″ = f(t, Derivative(2))

- alternatively, if one knows the index `i` associated to the segment of
  interest, then one can do

      s⃗  = f(i, ζ)
      s⃗′ = f(i, ζ, Derivative(1))
      s⃗″ = f(i, ζ, Derivative(2))

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

    One should use for example [`UnitTangent`](@ref) or [`CurvatureVector`](@ref) if one
    wants the derivatives with respect to the arc length ``ξ`` (which are more geometrically
    meaningful, and guaranteed to be orthogonal to each other).
    There is also [`normalise_derivatives`](@ref) which can be more efficient when one
    already has the derivatives with respect to ``t``.
"""
abstract type AbstractFilament{T} <: AbstractVector{Vec3{T}} end

Base.IndexStyle(::Type{<:AbstractFilament}) = IndexLinear()

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
    discretisation_method(f::AbstractFilament) -> DiscretisationMethod

Return the method used to discretise the filament based on its node locations.
"""
function discretisation_method end

"""
    nodes(f::AbstractFilament{T}) -> AbstractVector{T}

Return the nodes (or discretisation points) ``\\bm{s}_i`` of the filament.
"""
nodes(f::AbstractFilament) = f.Xs

@doc raw"""
    knots(f::AbstractFilament{T}) -> AbstractVector{T}

Return parametrisation knots ``t_i`` of the filament.

Filaments are parametrised by ``\bm{s}(t)`` for ``t ∈ [0, T]``.
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

Base.checkbounds(::Type{Bool}, f::AbstractFilament, I...) = checkbounds(Bool, nodes(f), I...)

"""
    Base.setindex!(f::AbstractFilament{T}, v, i::Int) -> Vec3{T}

Set coordinates of discretisation point ``\\bm{s}_i``.
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
include("closed_filament.jl")
include("segments.jl")
include("integrate.jl")
include("min_distance.jl")
include("refinement.jl")
include("reconnections.jl")
include("from_vector_field.jl")

include("interpolation/interpolation.jl")
include("interpolation/hermite.jl")

include("local/finitediff.jl")
include("local/closed_filament.jl")

include("spline/spline.jl")
include("spline/closed_filament.jl")

include("quantities.jl")
include("utils.jl")

include("plotting.jl")

"""
    Base.getindex(f::AbstractFilament{T}, i::Int) -> Vec3{T}
    Base.getindex(f::AbstractFilament{T}, i::Int, ::Derivative{n}) -> Vec3{T}
    Base.getindex(f::AbstractFilament{T}, i::Int, ::GeometricQuantity)

Return coordinates of discretisation point ``\\bm{s}_i``.

One may also obtain derivatives and other geometric quantities at point ``\\bm{s}_i``
by passing an optional [`Derivative`](@ref) or [`GeometricQuantity`](@ref).
"""
Base.@propagate_inbounds Base.getindex(f::AbstractFilament, i::Int) = nodes(f)[i]

Base.@propagate_inbounds Base.getindex(
    f::AbstractFilament, i::Int, d::Union{Derivative, GeometricQuantity},
) = f(AtNode(i), d)

# This is the default parametrisation for any generic discretisation method.
# Specific discretisation methods (e.g. Fourier) may choose a different default
# parametrisation.
default_parametrisation(::DiscretisationMethod) = (Xs, i) -> @inbounds(Xs[i])

@doc raw"""
    Filaments.init(
        ClosedFilament{T}, N::Integer, method::DiscretisationMethod;
        offset = zero(Vec3{T}),
        parametrisation = Filaments.default_parametrisation(method),
    ) -> ClosedFilament{T}

Allocate data for a closed filament with `N` discretisation points.

The element type `T` can be omitted, in which case the default `T = Float64` is used.

The optional `offset` keyword argument allows to define a filament with a spatial
offset between points `f[i]` and `f[i + N]`. By default the offset is zero,
meaning that the filament is a closed loop. This can be used for defining
infinite (so not really closed) filaments living in periodic domains.

Possible discretisation methods include [`CubicSplineMethod`](@ref) and
[`FiniteDiffMethod`](@ref).

# Customising the filament parametrisation

By default the filament is parametrised as ``\bm{s}(t)`` such that the parameter ``t``
roughly corresponds to the arc length ``ξ``.
More precisely, the discrete values ``t_i`` at the discretisation points are defined from
the distances between discretisation points:

```math
t_{i + 1} = t_{i} + |\bm{s}_{i + 1} - \bm{s}_i|, \quad t_1 = 0.
```

One can use the `parametrisation` keyword argument to change this.
This must be a function `parametrisation(Xs, i)` which takes the vector `Xs` of
discretisation points and the index `i` of the current point.
The default parametrisation function is simply `parametrisation(Xs, i) = Xs[i]`, which
corresponds to the definition above.

For instance, to parametrise filaments according to the ``z`` increments between
discretisation points, one can pass `parametrisation(Xs, i) = Xs[i].z`.

!!! warning

    Changing the filament parametrisation is **experimental**.
    This feature may change/disappear in the future.

"""
function init end

"""
    Filaments.init(
        ClosedFilament, points::AbstractVector{<:Vec3}, method::DiscretisationMethod;
        [kws...],
    )

Initialise new filament with the chosen discretisation points.

Note that [`update_coefficients!`](@ref) does not need to be called after using
this variant (until node locations change, of course!).
"""
function init(
        ::Type{FilamentType}, positions::AbstractVector{<:Vec3}, method::DiscretisationMethod;
        kws...,
    ) where {FilamentType <: AbstractFilament}
    @assert !(positions isa PaddedVector)  # avoid recursion...
    M = npad(method)
    N = length(positions)
    Xs = PaddedVector{M}(similar(positions, N + 2M))
    copyto!(Xs, positions)
    init(FilamentType, Xs, method; kws...)  # calls variant which takes a PaddedVector
end

"""
    update_coefficients!(f::AbstractFilament; knots = nothing)

Compute coefficients needed to perform inter-node interpolations and estimate
derivatives.

Uses the current locations of the filament nodes. If nodes change, this
function should be called to update the coefficients.

By default, this function also updates the parametrisation knots ``t_i``
according to the current node positions. One can override this by passing a
`knots` vector as a keyword argument. In particular, one can pass `knots = knots(f)` to
keep the parametrisation knots unchanged.

This function will fail if the number of filament nodes is smaller than that
required by the discretisation method. For instance, closed filaments
discretised using cubic splines must have at least 3 nodes.
"""
function update_coefficients! end

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
    normalise_derivatives(ṡ::Vec3, s̈::Vec3) -> (s⃗′, s⃗″)
    normalise_derivatives((ṡ, s̈)::NTuple)   -> (s⃗′, s⃗″)

Return derivatives with respect to the arc length ``ξ``, from derivatives with
respect to the parameter ``t``.

The returned derivatives satisfy:

- ``\\bm{s}' ≡ t̂`` is the **unit tangent** vector;

- ``\\bm{s}'' ≡ ρ n̂`` is the **curvature** vector, where ``n̂`` is the normal unit
  vector (with ``t̂ ⋅ n̂ = 0``) and ``ρ = R^{-1}`` is the curvature (and R the
  curvature radius).
"""
function normalise_derivatives(Ẋ::Vec3, Ẍ::Vec3)
    s′ = normalize(Ẋ)  # unit tangent vector
    s″ = (Ẍ - (Ẍ ⋅ s′) * s′) ./ sum(abs2, Ẋ)  # curvature vector
    s′, s″
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
function _update_knots_periodic!(
        parametrisation::F, ts::PaddedVector, Xs::PaddedVector, ::Nothing = nothing,
    ) where {F <: Function}
    @assert eachindex(ts) == eachindex(Xs)
    ts[begin] = 0
    inds = eachindex(ts)
    @assert npad(ts) == npad(Xs) ≥ 1
    yi = parametrisation(Xs, first(inds))
    @inbounds for i ∈ inds
        yj = parametrisation(Xs, i + 1)
        ts[i + 1] = ts[i] + norm(yj - yi)
        yi = yj
    end
    L = ts[end + 1] - ts[begin]  # knot period
    pad_periodic!(ts, L)
end

# In this case, override the computation of knots and copy the input knot vector.
function _update_knots_periodic!(
        ::F, ts::PaddedVector, Xs::PaddedVector, ts_in::PaddedVector,
    ) where {F}
    if ts !== ts_in  # avoid copy if input and output vectors are the same
        copyto!(ts, ts_in)  # fails if arrays have different dimensions
    end
    ts
end

# Similar to above but for a non-PaddedVector, so we need to apply periodic padding on the
# output.
function _update_knots_periodic!(
        ::F, ts::PaddedVector, Xs::PaddedVector, ts_in::AbstractVector,
    ) where {F}
    @assert !(ts_in isa PaddedVector)
    length(ts_in) == length(ts) + 1 ||
        throw(DimensionMismatch("input knots should include the endpoint"))
    L = ts_in[end] - ts[begin]  # knot period
    copyto!(ts, @view(ts_in[begin:end - 1]))
    pad_periodic!(ts, L)
    ts
end

"""
    redistribute_nodes!(f::AbstractFilament) -> f

Redistribute nodes of the filament so that they are (approximately) equally spaced.

More precisely, this function repositions the filament nodes such that the knot spacing
``t_{i + 1} - t_{i}`` is constant.
In other words, the new locations satisfy `f[i] = f((i - 1) * Δt)` where ``Δt = t_{N + 1} / N`` is
the knot spacing, ``N`` is the number of nodes, and the index ``N + 1`` refers to the
filament endpoint (which is equal to the starting point for a closed filament).
"""
function redistribute_nodes!(f::AbstractFilament)
    ta, tb = knotlims(f)
    N = length(f)
    ts = range(ta, tb; length = N + 1)  # should include the endpoint for `update_coefficients!`
    @assert eachindex(f) == eachindex(ts)[Base.OneTo(N)]
    @inbounds for i ∈ eachindex(f)
        # We interpolate at ts[i] from the old interpolation coefficients to get the new
        # node location f[i].
        # This assumes nodes and interpolation coefficients are stored in different places.
        f[i] = f(ts[i])
    end
    update_coefficients!(f; knots = ts)
    f
end

end
