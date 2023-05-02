using ..BasicTypes: Infinity

"""
    fold_periodic!(Xs::AbstractVector{<:Vec3}, Ls::NTuple{3, Real})

Fold a set of coordinates onto the main unit cell.

The idea is that the filament coordinates are (mainly) in the main unit cell, given by
``[0, L₁] × [0, L₂] × [0, L₃]``, where ``Lᵢ`` is the period in each direction.

This is nice for visualisations, and might also improve performance of the
folding required by long-range computations.

To avoid creating discontinuities in the filament, what this function actually
does is to make sure that the *average* among all filament nodes is in the main
unit cell.
"""
function fold_periodic!(Xs::AbstractVector{<:Vec3}, periods::NTuple{3, Real})
    Xmean = sum(Xs) ./ length(Xs)
    noffsets = map(_count_periodic_offsets, Xmean, periods)
    if all(==(0), noffsets)  # Xmean is already in the main unit cell
        return Xs
    end
    δx⃗ = oftype(Xmean, noffsets .* periods)
    for (i, x⃗) ∈ pairs(Xs)
        Xs[i] = x⃗ + δx⃗
    end
    Xs
end

"""
    fold_periodic!(f::AbstractFilament, Ls)

Fold filament nodes onto the main unit cell.
"""
function fold_periodic!(f::AbstractFilament, periods::NTuple{3, Real})
    fold_periodic!(nodes(f), periods)
    pad_periodic!(nodes(f))
    # TODO also update coefficients?
    f
end

function _count_periodic_offsets(x::Real, L::Real)
    offset = 0
    while x < 0
        x += L
        offset += 1
    end
    while x ≥ L
        x -= L
        offset -= 1
    end
    offset
end

_count_periodic_offsets(x::Real, ::Infinity) = 0

_find_knot_segment(ileft::Integer, tlims, ts, t) = (ileft, t)

function _find_knot_segment(::Nothing, tlims, ts, t)
    ta, tb = tlims
    T = tb - ta
    # This enables evaluation outside of the knot limits.
    while t ≥ tb
        t -= T
    end
    while t < ta
        t += T
    end
    i = searchsortedlast(ts, t) :: Int
    i, t
end
