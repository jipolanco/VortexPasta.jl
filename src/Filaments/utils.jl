using ..BasicTypes: Infinity

"""
    fold_periodic!(Xs::AbstractVector{<:Vec3}, Ls::NTuple{3, Real}) -> Bool

Fold a set of coordinates onto the main unit cell.

The idea is that the filament coordinates are (mainly) in the main unit cell, given by
``[0, L₁] × [0, L₂] × [0, L₃]``, where ``Lᵢ`` is the period in each direction.

This is nice for visualisations, and might also improve performance of the
folding required by long-range computations.

To avoid creating discontinuities in the filament, what this function actually
does is to make sure that the *average* among all filament nodes is in the main
unit cell.

Returns `true` if coordinates were modified, `false` otherwise.
In the first case, one may want to call [`update_coefficients!`](@ref) to update the curve
representation.
"""
function fold_periodic!(Xs::AbstractVector{<:Vec3}, periods::NTuple{3, Real})
    Xmean = sum(Xs) ./ length(Xs)
    noffsets = map(_count_periodic_offsets, Xmean, periods)
    if all(==(0), noffsets)  # Xmean is already in the main unit cell
        return false
    end
    δx⃗ = oftype(Xmean, noffsets .* periods)
    for (i, x⃗) ∈ pairs(Xs)
        Xs[i] = x⃗ + δx⃗
    end
    true
end

"""
    fold_periodic!(f::AbstractFilament, Ls) -> Bool

Fold filament nodes onto the main unit cell.

Returns `true` if coordinates were modified, `false` otherwise.

If `true`, coefficients may need to be updated afterwards using [`update_coefficients!`](@ref).
"""
function fold_periodic!(f::AbstractFilament, periods::NTuple{3, Real})
    modified = fold_periodic!(nodes(f), periods) :: Bool
    if modified
        pad_periodic!(nodes(f))
    end
    modified
end

fold_periodic!(f, periods::Vec3) = fold_periodic!(f, Tuple(periods))

# ======================================================================================== #

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

# ======================================================================================== #

# Return coordinate corresponding to x⃗ but in the main periodic cell.
# Here the periods are Ls = (Lx, Ly, Lz). If a component is `nothing`, then we do nothing in
# that direction (as if it was a non-periodic dimension).
to_main_periodic_cell(x⃗, Ls::Tuple) = oftype(x⃗, map(to_main_periodic_cell, x⃗, Ls))
to_main_periodic_cell(x, L::Nothing) = x  # don't do anything (non-periodic case)

function to_main_periodic_cell(x::T, L::Real) where {T}
    while x ≥ L
        x -= L
    end
    while x < 0
        x += L
    end
    x::T
end

# ======================================================================================== #

# Returns `true` if the distance between two points is larger than half the domain period in
# at least one direction. This is generally used to see if there is a periodic "jump"
# between two consecutive points of a filament (used for plots / VTK exports only).
is_jump(x⃗, x⃗_prev, Lhs::Tuple) = any(splat(is_jump), zip(x⃗, x⃗_prev, Lhs))
is_jump(x, x_prev, Lh::Nothing) = false
is_jump(x, x_prev, Lh::Real) = abs(x - x_prev) ≥ Lh

# ======================================================================================== #

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

# ======================================================================================== #

# TESTING / EXPERIMENTAL
# This function may be removed in the future.
# The idea is to update the parametrisation of `f` to follow more closely the
# actual arc lengths of the filament. Not sure if it's worth it...
function recompute_parametrisation!(f::ClosedFilament, quad::AbstractQuadrature)
    m = interpolation_method(f)
    _recompute_parametrisation!(m, f, quad)
end

# In the case of straight segments (linear interpolation), the parametrisation
# cannot be improved from its initial estimation.
_recompute_parametrisation!(::HermiteInterpolation{0}, f::AbstractFilament, quad) = f

function _recompute_parametrisation!(::Any, f::AbstractFilament, quad)
    (; ts,) = f
    @assert npad(ts) ≥ 1
    tnext = ts[begin]
    for i ∈ eachindex(ts)
        # Estimate arc length from ts[i] to ts[i + 1]
        ℓ = integrate(f, i, quad) do f, i, ζ
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
