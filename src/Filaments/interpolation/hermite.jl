export
    HermiteInterpolation

"""
    HermiteInterpolation{M} <: LocalInterpolationMethod

Hermite interpolation of continuity ``C^M`` at interpolation points.

Hermite interpolations are obtained using curve derivatives up to order ``M``.

## Allowed cases

- for ``M = 0`` this is simply linear interpolation (note that curvatures
  cannot be estimated from linear interpolations);

- for ``M = 1`` this is the standard Hermite interpolation (piecewise cubic
  polynomials, requires first derivatives at interpolation points);

- for ``M = 2`` this is a quintic Hermite interpolation requiring first and
  second derivatives at interpolation points.
"""
struct HermiteInterpolation{M} <: LocalInterpolationMethod end

HermiteInterpolation(M::Int) = HermiteInterpolation{M}()

continuity(::Type{HermiteInterpolation{M}}) where {M} = M
required_derivatives(::Type{HermiteInterpolation{M}}) where {M} = M
required_derivatives(m::HermiteInterpolation) = required_derivatives(typeof(m))

# Linear interpolation
# Coordinates and derivatives are expected to be normalised so that t ∈ [0, 1]
function interpolate(
        ::HermiteInterpolation{0}, ::Derivative{N},
        t::Number, Xs::NTuple{2}, etc...,
    ) where {N}
    N::Int
    if N === 0
        (1 - t) * Xs[1] + t * Xs[2]
    elseif N === 1
        Xs[2] - Xs[1]
    elseif N ≥ 2
        zero(Xs[1])
    end
end

# Cubic Hermite interpolation
@inline function interpolate(
        ::HermiteInterpolation{1}, ::Derivative{0},
        t::Number, Xs::NTuple{2}, Xs′::NTuple{2}, etc...,
    )
    t2 = t * t
    t3 = t2 * t
    (
        (2 * t3 - 3 * t2 + 1) * Xs[1]
        +
        (-2 * t3 + 3 * t2) * Xs[2]
        +
        (t3 - 2 * t2 + t) * Xs′[1]
        +
        (t3 - t2) * Xs′[2]
    )
end

@inline function interpolate(
        ::HermiteInterpolation{1}, ::Derivative{1},
        t::Number, Xs::NTuple{2}, Xs′::NTuple{2}, etc...,
    )
    t2 = t * t
    (
        (6 * t2 - 6 * t) * Xs[1]
        +
        (-6 * t2 + 6 * t) * Xs[2]
        +
        (3 * t2 - 4 * t + 1) * Xs′[1]
        +
        (3 * t2 - 2t) * Xs′[2]
    )
end

@inline function interpolate(
        ::HermiteInterpolation{1}, ::Derivative{2},
        t::Number, Xs::NTuple{2}, Xs′::NTuple{2}, etc...,
    )
    (
        (12 * t - 6) * Xs[1]
        +
        (-12 * t + 6) * Xs[2]
        +
        (6 * t - 4) * Xs′[1]
        +
        (6 * t - 2) * Xs′[2]
    )
end

# Quintic Hermite interpolation
@inline function interpolate(
        ::HermiteInterpolation{2}, ::Derivative{0},
        t::Number, Xs::NTuple{2}, Xs′::NTuple{2}, Xs″::NTuple{2}, etc...,
    )
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    t5 = t2 * t3
    (
        (1 - 10t3 + 15t4 - 6t5) * Xs[1]
        +
        (10t3 - 15t4 + 6t5) * Xs[2]
        +
        (t - 6t3 + 8t4 - 3t5) * Xs′[1]
        +
        (-4t3 + 7t4 - 3t5) * Xs′[2]
        +
        (t2 - 3t3 + 3t4 - t5) / 2 * Xs″[1]
        +
        (t3 - 2t4 + t5) / 2 * Xs″[2]
    )
end

@inline function interpolate(
        ::HermiteInterpolation{2}, ::Derivative{1},
        t::Number, Xs::NTuple{2}, Xs′::NTuple{2}, Xs″::NTuple{2}, etc...,
    )
    t2 = t * t
    t3 = t2 * t
    t4 = t2 * t2
    (
        30 * (-t2 + 2t3 - t4) * Xs[1]
        +
        30 * (t2 - 2t3 + t4) * Xs[2]
        +
        (1 - 18t2 + 32t3 - 15t4) * Xs′[1]
        +
        (-12t2 + 28t3 - 15t4) * Xs′[2]
        +
        (2t - 9t2 + 12t3 - 5t4) / 2 * Xs″[1]
        +
        (3t2 - 8t3 + 5t4) / 2 * Xs″[2]
    )
end

@inline function interpolate(
        ::HermiteInterpolation{2}, ::Derivative{2},
        t::Number, Xs::NTuple{2}, Xs′::NTuple{2}, Xs″::NTuple{2}, etc...,
    )
    t2 = t * t
    t3 = t2 * t
    (
        (-60t + 180t2 - 120t3) * Xs[1]
        +
        (60t - 180t2 + 120t3) * Xs[2]
        +
        (-36t + 96t2 - 60t3) * Xs′[1]
        +
        (-24t + 84t2 - 60t3) * Xs′[2]
        +
        (2 - 18t + 36t2 - 20t3) / 2 * Xs″[1]
        +
        (6t - 24t2 + 20t3) / 2 * Xs″[2]
    )
end

function interpolate(method::HermiteInterpolation, d::Derivative, args...)
    throw(ArgumentError(
        lazy"interpolation of $d with $method not currently implemented"
    ))
end

## ================================================================================ ##

# Evaluate interpolation using local parametrisation ζ ∈ [0, 1] (within segment `i`)
function evaluate(
        method::HermiteInterpolation{M}, coefs::DiscretisationCoefs, ts::PaddedVector,
        i::Int, ζ::Number, deriv::Derivative{N} = Derivative(0),
    ) where {M, N}
    (; cs, cderivs,) = coefs
    checkbounds(ts, i)
    @assert npad(cs) ≥ 1
    @assert all(X -> npad(X) ≥ 1, cderivs)
    @inbounds ℓ_i = ts[i + 1] - ts[i]
    @assert M ≤ 2
    α = 1 / (ℓ_i^N)  # normalise returned derivative by ℓ^N
    values_full = (cs, cderivs...)
    ℓ_norm = Ref(one(ℓ_i))
    values_i = ntuple(Val(M + 1)) do m
        @inline
        X = values_full[m]
        @inbounds data = (X[i], X[i + 1]) .* ℓ_norm  # normalise m-th derivative by ℓ^m
        ℓ_norm[] *= ℓ_i
        data
    end
    α * interpolate(method, deriv, ζ, values_i...)
end

# Evaluate interpolation using global parametrisation t ∈ [0, T]
function evaluate(
        method::HermiteInterpolation, coefs::DiscretisationCoefs, ts::PaddedVector,
        t_in::Number, deriv::Derivative = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    )
    ta = ts[begin]
    tb = ts[end + 1]
    i, t = _find_knot_segment(ileft, (ta, tb), ts, t_in)
    ζ = @inbounds (t - ts[i]) / (ts[i + 1] - ts[i])
    evaluate(method, coefs, ts, i, ζ, deriv)
end

## ================================================================================ ##
# Interpolation of filament coordinates and derivatives

function _interpolate(
        m::HermiteInterpolation, f::ClosedFilament, t_in::Number, d::Derivative = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    )
    (; ts,) = f
    i, t = _find_knot_segment(ileft, knotlims(f), ts, t_in)
    ζ = @inbounds (t - ts[i]) / (ts[i + 1] - ts[i])
    y = evaluate(m, f.coefs, ts, i, ζ, d)
    _deperiodise_hermite(d, y, f, t, t_in)
end

function _interpolate(
        method::HermiteInterpolation, f::ClosedFilament,
        i::Int, ζ::Number, deriv::Derivative,
    )
    evaluate(method, f.coefs, f.ts, i, ζ, deriv)
end

function _deperiodise_hermite(::Derivative{0}, y, f::ClosedFilament, t, t_in)
    (; Xoffset,) = f
    Xoffset === zero(Xoffset) && return y
    t == t_in && return y
    ta, tb = knotlims(f)
    T = tb - ta
    dt = t_in - t  # expected to be a multiple of T
    y + dt / T * Xoffset
end

_deperiodise_hermite(::Derivative, y, args...) = y  # derivatives (n ≥ 1): shift not needed
