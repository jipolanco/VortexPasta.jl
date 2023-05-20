using LinearAlgebra: norm
using StaticArrays: SVector

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
    NoRefinement <: RefinementCriterion
    NoRefinement()

Used to disable filament refinement.
"""
struct NoRefinement <: RefinementCriterion end

@inline Base.getproperty(crit::NoRefinement, name::Symbol) = _getproperty(crit, Val(name))
_getproperty(::NoRefinement, ::Val{:inds}) = SVector{0, Int}()
_getproperty(::NoRefinement, ::Val{:remove}) = SVector{0, Bool}()

_nodes_to_refine!(::AbstractFilament, ::NoRefinement) = (0, 0)

"""
    BasedOnCurvature <: RefinementCriterion
    BasedOnCurvature(ρℓ_max::Real, ρℓ_min = ρℓ_max / 2.5; ℓ_max = Inf)

Curvature-based refinement criterion.

## Node insertion

According to this criterion, a filament is locally refined if:

```math
ρ \\left|\\bm{X}_{i + 1} - \\bm{X}_{i}\\right| > (ρℓ)_{\\mathrm{max}}
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
To limit this, set the keyword argument `ℓ_max` to some finite value
determining the maximum length of a segment.
This setting will be roughly respected when removing nodes.
"""
struct BasedOnCurvature <: RefinementCriterion
    ρℓ_max :: Float64
    ρℓ_min :: Float64
    ℓ_max  :: Float64
    inds   :: Vector{Int}   # indices of nodes or segments to modify
    remove :: Vector{Bool}  # determine whether nodes should be removed or segments should be refined
    BasedOnCurvature(ρℓ_max, ρℓ_min; ℓ_max = Inf) = new(ρℓ_max, ρℓ_min, ℓ_max, Int[], Bool[])
end

BasedOnCurvature(ρℓ_max; kws...) = BasedOnCurvature(ρℓ_max, ρℓ_max / 2.5; kws...)

function _nodes_to_refine!(f::AbstractFilament, crit::BasedOnCurvature)
    (; ρℓ_max, ρℓ_min, ℓ_max, inds, remove,) = crit
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
        # ℓ = ts[i + 1] - ts[i]  # assumes parametrisation corresponds to node distance
        ℓ = norm(f[i + 1] - f[i])
        ρ = (f[i, CurvatureScalar()] + f[i + 1, CurvatureScalar()]) / 2
        # ρ_alt = f(i, 0.5, CurvatureScalar())  # this is likely more expensive, and less accurate for FiniteDiff
        ρℓ = ρ * ℓ
        if ρℓ > ρℓ_max
            push!(inds, i)
            push!(remove, false)
            n_add += 1
        elseif ρℓ < ρℓ_min && ℓ < ℓ_max / 2  # so that the new ℓ is roughly smaller than ℓ_max
            push!(inds, i + 1)
            push!(remove, true)
            n_rem += 1
            skipnext = true
        end
    end
    n_add, n_rem
end
