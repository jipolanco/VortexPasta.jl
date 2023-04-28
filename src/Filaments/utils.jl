"""
    fold_periodic!(Xs::AbstractVector{<:Vec3}, Ls::NTuple{3, Real})

Fold a set of coordinates onto the main unit cell.

The idea is that the filament coordinates are (mainly) in the main unit cell, given by
``[0, L₁] × [0, L₂] × [0, L₃]``, where ``Lᵢ`` is the period in each direction.

This is nice for visualisations, and might also improve performance of the
folding required by long-range computations (see [`BiotSavart.folding_limits`](@ref)).

To avoid creating discontinuities in the filament, what this function actually
does is to make sure that the *average* among all filament nodes is in the main
unit cell.
"""
function fold_periodic!(Xs::AbstractVector{<:Vec3}, periods::NTuple{3, Real})
    Xmean = sum(Xs) ./ length(Xs)
    noffsets = map(Xmean, periods) do x, L
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
    if all(==(0), noffsets)  # Xmean is already in the main unit cell
        return Xs
    end
    δx⃗ = oftype(Xmean, noffsets .* periods)
    for (i, x⃗) ∈ pairs(Xs)
        Xs[i] = x⃗ + δx⃗
    end
    Xs
end
