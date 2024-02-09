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
