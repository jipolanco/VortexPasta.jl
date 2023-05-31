using ..BasicTypes: Zero

"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`BasedOnDistance`](@ref): reconnects filament segments which are closer than a critical distance.
"""
abstract type ReconnectionCriterion end

"""
    distance(crit::ReconnectionCriterion) -> Real

Return the critical distance associated to the reconnection criterion.
"""
function distance end

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = Zero()

"""
    BasedOnDistance <: ReconnectionCriterion
    BasedOnDistance(d_crit; cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

The keyword argument `cos_max` disables the reconnection of nearly parallel segments.
The default value `cos_max = 0.97` disables reconnections when the angle between lines
is ``θ < \\arccos(0.97) ≈ 14°``.
"""
struct BasedOnDistance <: ReconnectionCriterion
    dist    :: Float64
    cos_max :: Float64

    function BasedOnDistance(dist; cos_max = 0.97)
        new(dist, cos_max)
    end
end

distance(c::BasedOnDistance) = c.dist

"""
    split!(f::ClosedFilament, i::Int, j::Int) -> (f₁, f₂)

Split closed filament into two filaments.

Assuming `j > i`, the resulting filaments are composed of nodes `f[i + 1:j]` and
`f[(j + 1:end) ∪ (begin:i)]`.

In practice, a split makes sense when the nodes `f[i]` and `f[j]` are "close".

Note that the filament `f` is modified by this function, and is returned as the filament `f₁`.

One should generally call [`update_coefficients!`](@ref) on both filaments after a split.
"""
function split!(f::ClosedFilament, i::Int, j::Int)
    i ≤ j || return split!(f, j, i)

    n1 = j - i
    n2 = length(f) - n1

    f1 = f
    f2 = similar(f, n2)

    # Fill f2
    l = firstindex(f2) - 1
    for k ∈ (j + 1):lastindex(f)
        f2[l += 1] = f[k]
    end
    for k ∈ firstindex(f):i
        f2[l += 1] = f[k]
    end
    @assert l == n2

    # Shift values and then resize f1
    l = firstindex(f1) - 1
    for k ∈ (i + 1):j
        f1[l += 1] = f[k]
    end
    @assert l == n1
    resize!(f1, n1)

    f1, f2
end
