export NoRefinement,
       RefineBasedOnCurvature,
       RefineBasedOnSegmentLength

using LinearAlgebra: norm
using StaticArrays: SVector

"""
    insert_node!(f::AbstractFilament, i::Integer, [ζ = 0.5]) -> Vec3

Insert node in-between locations `f[i]` and `f[i + 1]`.

The optional argument `ζ ∈ [0, 1]` corresponds to the relative location of the new node
within the segment. By default it is set to `ζ = 0.5`, which corresponds to an estimation of
the middle of the segment.

Note that [`update_after_changing_nodes!`](@ref) must be called after inserting one or more
nodes.

See also [`remove_node!`](@ref).

Returns the inserted node.
"""
function insert_node! end

insert_node!(f::AbstractFilament, i::Integer) = insert_node!(f, i, 0.5)

"""
    remove_node!(f::AbstractFilament, i::Integer) -> Vec3

Remove node at location `f[i]`.

Note that [`update_after_changing_nodes!`](@ref) must be called after removing one or more
nodes.

See also [`insert_node!`](@ref).

Returns the removed node.
"""
function remove_node! end

"""
    update_after_changing_nodes!(f::AbstractFilament)

Update filament fields after changing nodes.

Depending on the filament discretisation method, this can recompute derivatives, knots or
discretisation coefficients.

Should be called after inserting or removing filament nodes.
See [`insert_node!`](@ref) and [`remove_node!`](@ref).
"""
function update_after_changing_nodes! end

"""
    RefinementCriterion

Abstract type describing a curve refinement criterion.

Implemented refinement criteria are:

- [`NoRefinement`](@ref): disables refinement;

- [`RefineBasedOnSegmentLength`](@ref): enforces a minimum and maximum distance between
  neighbouring filament nodes;

- [`RefineBasedOnCurvature`](@ref): inserts more nodes on highly-curved filament segments.

"""
abstract type RefinementCriterion end

function _nodes_to_refine!(f::AbstractFilament, crit::RefinementCriterion)
    (; cache,) = crit
    (; inds_add, inds_rem,) = cache
    ts = knots(f)
    empty!(cache)
    skipnext = false
    iter = eachindex(segments(f))  # iterate over segments of the unmodified filament
    for i ∈ iter
        if skipnext
            skipnext = false
            continue
        end
        action = _refinement_action(crit, f, i)
        if action === :insert
            @debug lazy"Inserting node at t = $((ts[i] + ts[i + 1]) / 2)"
            push!(inds_add, i)
        elseif action === :remove
            @debug lazy"Removing node at t = $(ts[i + 1])"
            push!(inds_rem, i + 1)
            skipnext = true
        end
    end
    cache
end

"""
    refine!(f::AbstractFilament, crit::RefinementCriterion) -> (Int, Int)

Refine the filament according to a given criterion.

More precisely, this function can add and remove discretisation points according to the
chosen refinement criterion.

Returns the number of added and removed nodes.

Example usage:

    crit = RefineBasedOnCurvature(0.5)
    refine!(f, crit)

"""
function refine!(f::AbstractFilament, crit::RefinementCriterion)
    method = discretisation_method(f)
    _refine!(method, f, crit)
end

# Default implementation
function _refine!(method::DiscretisationMethod, f, crit)
    @assert method === discretisation_method(f)
    N = length(f)  # original number of nodes

    # Determine where to add or remove nodes.
    cache = _nodes_to_refine!(f, crit)
    (; inds_add, inds_rem,) = cache
    n_add = length(inds_add)
    n_rem = length(inds_rem)
    iszero(n_add + n_rem) && return (n_add, n_rem)

    # Note: we add all knots first, then we remove all knots to be removed.
    if n_add > 0
        sizehint!(f, N + n_add)
    end

    # 1. Add nodes.
    # We iterate in reverse to avoiding the need to shift indices (assuming indices are
    # sorted).
    for i ∈ reverse(inds_add)
        insert_node!(f, i)
        # Shift indices to be removed if needed
        for (n, j) ∈ pairs(inds_rem)
            if j > i
                inds_rem[n] += 1
            end
        end
    end

    # 2. Remove nodes
    for i ∈ reverse(inds_rem)
        remove_node!(f, i)
    end

    if n_add + n_rem > 0
        @assert length(nodes(f)) == N + n_add - n_rem
        update_after_changing_nodes!(f; removed = n_rem > 0)
    end

    n_add, n_rem
end

struct RefinementCache
    inds_rem :: Vector{Int}  # indices of nodes to remove
    inds_add :: Vector{Int}  # indices of segments to refine
    RefinementCache() = new(Int[], Int[])
end

function Base.empty!(c::RefinementCache)
    empty!(c.inds_rem)
    empty!(c.inds_add)
    c
end

"""
    NoRefinement <: RefinementCriterion
    NoRefinement()

Used to disable filament refinement.
"""
struct NoRefinement <: RefinementCriterion end

@inline Base.getproperty(crit::NoRefinement, name::Symbol) = _getproperty(crit, Val(name))

# Return a fake cache.
function _getproperty(::NoRefinement, ::Val{:cache})
    inds_add = SVector{0, Int}()
    inds_rem = SVector{0, Int}()
    (; inds_add, inds_rem,)
end

_nodes_to_refine!(::AbstractFilament, crit::NoRefinement) = crit.cache

"""
    RefineBasedOnCurvature <: RefinementCriterion
    RefineBasedOnCurvature(ρℓ_max::Real, ρℓ_min = ρℓ_max / 2.5; ℓ_min = 0.0, ℓ_max = Inf)

Curvature-based refinement criterion.

## Node insertion

According to this criterion, a filament is locally refined if:

```math
ρ \\left|\\bm{X}_{i + 1} - \\bm{X}_{i}\\right| > (ρℓ)_{\\max}
```

where ``ρ`` is some estimation of the curvature of segment ``[i, i + 1]``.

For splines, this uses a classical knot insertion algorithm which preserves the
shape of the curve.

## Node removal

Similarly, filaments nodes are removed based on the value of `ρℓ_min`.
This value should be less than `ρℓ_max / 2` to avoid alternatively adding and
removing nodes when repeatedly calling [`refine!`](@ref).

For safety, two adjacent nodes will never be removed in a single call to `refine!`.

Note that, when filaments are nearly straight, this may lead to the removal of most nodes.
To limit this, set the keyword argument `ℓ_max` to some finite value determining the maximum
length of a segment.
Similarly, the keyword argument `ℓ_min` sets a lower limit for the distance between
neighbouring nodes.
"""
struct RefineBasedOnCurvature <: RefinementCriterion
    ρℓ_max :: Float64
    ρℓ_min :: Float64
    ℓ_max  :: Float64
    ℓ_min  :: Float64
    cache  :: RefinementCache
    RefineBasedOnCurvature(ρℓ_max, ρℓ_min; ℓ_max = Inf, ℓ_min = 0.0) =
        new(ρℓ_max, ρℓ_min, ℓ_max, ℓ_min, RefinementCache())
end

RefineBasedOnCurvature(ρℓ_max; kws...) = RefineBasedOnCurvature(ρℓ_max, ρℓ_max / 2.5; kws...)

function Base.show(io::IO, c::RefineBasedOnCurvature)
    (; ρℓ_max, ρℓ_min, ℓ_max, ℓ_min,) = c
    print(io, "RefineBasedOnCurvature($ρℓ_max, $ρℓ_min; ℓ_max = $ℓ_max, ℓ_min = $ℓ_min)")
end

function _refinement_action(crit::RefineBasedOnCurvature, f::AbstractFilament, i::Integer)
    (; ρℓ_min, ρℓ_max, ℓ_min, ℓ_max,) = crit
    ℓ = norm(f[i + 1] - f[i])
    ρ = (f[i, CurvatureScalar()] + f[i + 1, CurvatureScalar()]) / 2
    # ρ_alt = f(i, 0.5, CurvatureScalar())  # this is likely more expensive, and less accurate for FiniteDiff
    ρℓ = ρ * ℓ
    if ρℓ > ρℓ_max && ℓ > 2 * ℓ_min  # so that the new ℓ is roughly larger than ℓ_min
        :insert
    elseif ρℓ < ρℓ_min && ℓ < ℓ_max / 2  # so that the new ℓ is roughly smaller than ℓ_max
        :remove
    else
        :nothing
    end
end

"""
    RefineBasedOnSegmentLength <: RefinementCriterion
    RefineBasedOnSegmentLength(ℓ_min, ℓ_max = 2 * ℓ_min)

Refinement criterion imposing a minimum segment length.

This refinement criterion imposes neighbouring filament nodes to be at a distance
``ℓ ∈ [ℓ_{\\min}, ℓ_{\\max}]``. This means that:

- nodes are **inserted** if the distance between two nodes is ``ℓ > ℓ_{\\max}``.
  The insertion is done at an intermediate position using the functional representation of
  the filament (e.g. splines or Hermite interpolation);

- nodes are **removed** if the distance between two nodes is ``ℓ < ℓ_{\\min}``.
  For a filament which is strongly curved at that point, this means that local information
  is lost and that the filament is smoothed.
"""
struct RefineBasedOnSegmentLength <: RefinementCriterion
    ℓ_min :: Float64
    ℓ_max :: Float64
    cache :: RefinementCache
    function RefineBasedOnSegmentLength(ℓ_min, ℓ_max = 2 * ℓ_min)
        ℓ_min < ℓ_max || error(lazy"ℓ_min should be smaller than ℓ_max (got ℓ_max/ℓ_min = $ℓ_max/$ℓ_min)")
        new(ℓ_min, ℓ_max, RefinementCache())
    end
end

function Base.show(io::IO, c::RefineBasedOnSegmentLength)
    (; ℓ_min, ℓ_max,) = c
    print(io, "RefineBasedOnSegmentLength($ℓ_min, $ℓ_max)")
end

function _refinement_action(crit::RefineBasedOnSegmentLength, f::AbstractFilament, i::Integer)
    (; ℓ_min, ℓ_max,) = crit
    ℓ = norm(f[i + 1] - f[i])
    if ℓ > ℓ_max
        :insert
    elseif ℓ < ℓ_min
        :remove
    else
        :nothing
    end
end
