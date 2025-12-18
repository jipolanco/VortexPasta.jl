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

# Action to be performed for each filament node.
@enum RefinementAction::UInt8 begin
    REFINEMENT_DO_NOTHING = 0x00    # do nothing with the current node
    REFINEMENT_INSERT = 0x01 << 0   # insert a node after the current node
    REFINEMENT_REMOVE = 0x01 << 1   # remove the current node
end

function _nodes_to_refine!(
        f::AbstractFilament, crit::RefinementCriterion;
        actions = Memory{UInt8}(undef, length(f) + 1)  # includes endpoint
    )
    @assert eachindex(f) == eachindex(actions)[begin:(end - 1)]
    # (1) Compute actions to be performed
    for i in eachindex(f)
        actions[i] = Integer(_refinement_action(crit, f, i)::RefinementAction)
    end
    actions[end] = Integer(REFINEMENT_DO_NOTHING)  # endpoint
    # (2) Post-process actions:
    # - move removal operations one node to the right
    # - make sure we don't remove contiguous nodes
    # - count number of operations
    n_add = 0
    n_rem = 0
    removed_latest = false
    @inbounds for i in eachindex(f)
        action = actions[i]
        n_add += !iszero(action & Integer(REFINEMENT_INSERT))
        if removed_latest
            removed_latest = false
        elseif !iszero(action & Integer(REFINEMENT_REMOVE))
            actions[i] = action ⊻ Integer(REFINEMENT_REMOVE)  # unset REFINEMENT_REMOVE operation from node i
            actions[i + 1] = action | Integer(REFINEMENT_REMOVE)  # set REFINEMENT_REMOVE operation on node i + 1
            n_rem += 1
            removed_latest = true
        end
    end
    # Wrap things at the endpoint / startpoint (this may add a remove operation to the startpoint)
    actions[begin] = actions[begin] | actions[end]
    (; actions, n_add, n_rem)
end

"""
    refine!(f::AbstractFilament, crit::RefinementCriterion) -> (n_added, n_removed)

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

    buf = Bumper.default_buffer()

    @no_escape buf begin
        # Determine where to add or remove nodes.
        actions = @alloc(UInt8, length(f) + 1)
        (; n_add, n_rem) = _nodes_to_refine!(f, crit; actions)

        if !iszero(n_add + n_rem)
            # Temporary vector where new nodes will be stored
            Xs = nodes(f)
            ts = knots(f)
            V = eltype(Xs)
            T = eltype(ts)
            N = length(f)  # original number of nodes
            N_after = N + n_add - n_rem
            Xs_after = @alloc(V, N_after)
            ts_after = @alloc(T, N_after + 1)  # must include the endpoint

            i = firstindex(Xs) - 1  # current index in Xs, f, actions
            i_after = i             # current index in Xs_after
            @inbounds while i_after < lastindex(Xs_after)
                # @assert i < lastindex(Xs)
                i += 1
                i_after += 1
                action = actions[i]
                # Note: a node can have both an insert and a remove action, which corresponds to
                # action = 0x01 | 0x02 = 0x03. This means that the node itself will be removed, but
                # that a node will also be inserted on the segment to its right.
                insert = !iszero(action & Integer(REFINEMENT_INSERT))
                remove = !iszero(action & Integer(REFINEMENT_REMOVE))
                if insert
                    if !remove
                        # 1. Copy node i as usual
                        Xs_after[i_after] = Xs[i]
                        ts_after[i_after] = ts[i]
                        i_after += 1  # increment counter for second point
                    end
                    # 2. Insert second node in the middle of segment [i, i + 1] (or replace
                    # existent point if remove == true).
                    # Note: for splines we could use standard knot insertion algorithms
                    # (which are actually implemented), but that's probably not worth it,
                    # since in most cases we need to recompute coefficients later anyways.
                    X_new = f(i, T(0.5))  # interpolate position in the middle of the next segment
                    Xs_after[i_after] = X_new
                    ts_after[i_after] = T(0.5) * (ts[i] + ts[i + 1])
                elseif remove
                    # Skip copying node i, and decrement i_after since it will be copied in the next iteration.
                    i_after -= 1
                else  # if action == REFINEMENT_DO_NOTHING
                    # Simply copy node i
                    Xs_after[i_after] = Xs[i]
                    ts_after[i_after] = ts[i]
                end
            end

            ts_after[end] = ts[end + 1]  # endpoint

            # Now copy new nodes onto original filament
            resize!(f, N_after)
            for i in eachindex(f, Xs_after)
                @inbounds f[i] = Xs_after[i]
            end

            # Finally, recompute coefficients, preserving computed knots.
            if check_nodes(Bool, f)  # if we have enough points (otherwise filament should be deleted later)
                update_coefficients!(f; knots = ts_after)
            end
        end
    end

    n_add, n_rem
end

"""
    NoRefinement <: RefinementCriterion
    NoRefinement()

Used to disable filament refinement.
"""
struct NoRefinement <: RefinementCriterion end

refine!(::AbstractFilament, ::NoRefinement) = (0, 0)  # = (n_added, n_removed)

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
    RefineBasedOnCurvature(ρℓ_max, ρℓ_min; ℓ_max = Inf, ℓ_min = 0.0) =
        new(ρℓ_max, ρℓ_min, ℓ_max, ℓ_min)
end

RefineBasedOnCurvature(ρℓ_max; kws...) = RefineBasedOnCurvature(ρℓ_max, ρℓ_max / 2.5; kws...)

function Base.show(io::IO, c::RefineBasedOnCurvature)
    (; ρℓ_max, ρℓ_min, ℓ_max, ℓ_min,) = c
    print(io, "RefineBasedOnCurvature($ρℓ_max, $ρℓ_min; ℓ_max = $ℓ_max, ℓ_min = $ℓ_min)")
end

function _refinement_action(crit::RefineBasedOnCurvature, f::AbstractFilament, i::Integer)::RefinementAction
    (; ρℓ_min, ρℓ_max, ℓ_min, ℓ_max,) = crit
    ℓ = norm(f[i + 1] - f[i])
    ρ = (f[i, CurvatureScalar()] + f[i + 1, CurvatureScalar()]) / 2
    # ρ_alt = f(i, 0.5, CurvatureScalar())  # this is likely more expensive, and less accurate for FiniteDiff
    ρℓ = ρ * ℓ
    if ρℓ > ρℓ_max && ℓ > 2 * ℓ_min  # so that the new ℓ is roughly larger than ℓ_min
        REFINEMENT_INSERT
    elseif ρℓ < ρℓ_min && ℓ < ℓ_max / 2  # so that the new ℓ is roughly smaller than ℓ_max
        REFINEMENT_REMOVE
    else
        REFINEMENT_DO_NOTHING
    end
end

"""
    RefineBasedOnSegmentLength <: RefinementCriterion
    RefineBasedOnSegmentLength(ℓ_min, ℓ_max = 2 * ℓ_min; ρℓ_max = Inf)

Refinement criterion imposing a minimum segment length.

This refinement criterion imposes neighbouring filament nodes to be at a distance
``ℓ ∈ [ℓ_{\\min}, ℓ_{\\max}]``. This means that:

- nodes are **inserted** if the distance between two nodes is ``ℓ > ℓ_{\\max}``.
  The insertion is done at an intermediate position using the functional representation of
  the filament (e.g. splines or Hermite interpolation);

- nodes are **removed** if the distance between two nodes is ``ℓ < ℓ_{\\min}``.
  For a filament which is strongly curved at that point, this means that local information
  is lost and that the filament is smoothed.

The optional argument `ρℓ_max` sets the curvature limit so that higher curvature regions will be smoothed
and total line length is decreased:

- nodes are **removed** if the local normalised curvature is ``ρℓ > (ρℓ)_{\\max}``.
"""
struct RefineBasedOnSegmentLength <: RefinementCriterion
    ℓ_min   :: Float64
    ℓ_min²  :: Float64
    ℓ_max   :: Float64
    ℓ_max²  :: Float64
    ρℓ_max  :: Float64
    ρℓ_max² :: Float64
    function RefineBasedOnSegmentLength(ℓ_min, ℓ_max = 2 * ℓ_min; ρℓ_max = Inf)
        ℓ_min < ℓ_max || error(lazy"ℓ_min should be smaller than ℓ_max (got ℓ_max/ℓ_min = $ℓ_max/$ℓ_min)")
        new(ℓ_min, ℓ_min^2, ℓ_max, ℓ_max^2, ρℓ_max, ρℓ_max^2)
    end
end

function Base.show(io::IO, c::RefineBasedOnSegmentLength)
    (; ℓ_min, ℓ_max, ρℓ_max ) = c
    print(io, "RefineBasedOnSegmentLength($ℓ_min, $ℓ_max; ρℓ_max = $ρℓ_max)")
end

function _refinement_action(crit::RefineBasedOnSegmentLength, f::AbstractFilament, i::Integer)::RefinementAction
    (; ℓ_min², ℓ_max², ρℓ_max²) = crit
    ℓ² = sum(abs2, f[i + 1] - f[i])
    if ℓ² > ℓ_max²
        REFINEMENT_INSERT
    elseif ℓ² < ℓ_min²
        REFINEMENT_REMOVE
    elseif !isinf(ρℓ_max²) # if removing large curvature regions 
        ρ = (f[i, CurvatureScalar()] + f[i + 1, CurvatureScalar()]) / 2
        ρℓ² = ρ^2 * ℓ²
        if ρℓ² > ρℓ_max²
            REFINEMENT_REMOVE
        else 
            REFINEMENT_DO_NOTHING
        end
    else
        REFINEMENT_DO_NOTHING
    end
end
