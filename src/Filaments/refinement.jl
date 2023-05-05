"""
    RefinementCriterion

Abstract type describing a curve refinement criterion.
"""
abstract type RefinementCriterion end

"""
    refine!(f::AbstractFilament, crit::RefinementCriterion) -> (Int, Int)

Refine the filament according to a given criterion.

By default, it also removes nodes

Returns the number of added and removed nodes.

Example usage:

    crit = BasedOnCurvature(0.5)
    refine!(f, crit)

"""
function refine! end

"""
    BasedOnCurvature <: RefinementCriterion
    BasedOnCurvature(ρℓ_max::Real, ρℓ_min = ρℓ_max / 2)

Curvature-based refinement criterion.

According to this criterion, a filament is locally refined if:

```math
ρ_{i + 1/2} \\left|\\bm{X}_{i + 1} - \\bm{X}_{i}\\right| > (ρℓ)_{\\mathrm{max}}
```

where ``ρ_{i + 1/2}`` is the curvature in the middle of the segment.

Similarly, filaments nodes are removed based on the value of `ρℓ_min`.
"""
struct BasedOnCurvature <: RefinementCriterion
    ρℓ_max :: Float64
    ρℓ_min :: Float64
    inds   :: Vector{Int}   # indices of nodes or segments to modify
    remove :: Vector{Bool}  # determine whether nodes should be removed or segments should be refined
    BasedOnCurvature(ρℓ_max, ρℓ_min) = new(ρℓ_max, ρℓ_min, Int[], Bool[])
end

BasedOnCurvature(ρℓ_max) = BasedOnCurvature(ρℓ_max, ρℓ_max / 2.1)

function _nodes_to_refine!(f::AbstractFilament, crit::BasedOnCurvature)
    (; ρℓ_max, ρℓ_min, inds, remove,) = crit
    ts = knots(f)
    n_add = n_rem = 0
    empty!(inds)
    empty!(remove)
    skipnext = false
    iter = eachindex(segments(f))  # iterate over segments of the unmodified filament
    for i ∈ iter
        if skipnext
            skipnext = false
            continue
        end
        ℓ = ts[i + 1] - ts[i]  # assume parametrisation corresponds to node distance
        ρ = f(i, 0.5, CurvatureScalar())
        ρℓ = ρ * ℓ
        if ρℓ > ρℓ_max
            push!(inds, i)
            push!(remove, false)
            n_add += 1
        elseif ρℓ < ρℓ_min
            push!(inds, i + 1)
            push!(remove, true)
            n_rem += 1
            skipnext = true
        end
    end
    n_add, n_rem
end
