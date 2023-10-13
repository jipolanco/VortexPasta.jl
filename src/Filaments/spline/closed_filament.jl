function _update_coefficients_only!(
        ::CubicSplineMethod, f::ClosedFilament;
        only_derivatives = false,
    )
    (; ts, Xs, Xoffset,) = f
    (; cs, cderivs,) = f.coefs
    if !only_derivatives
        solve_cubic_spline_coefficients!(cs, ts, Xs; Xoffset,)
    end
    spline_derivative!(cderivs[1], cs, ts, Val(4))
    spline_derivative!(cderivs[2], cderivs[1], ts, Val(3))
    f
end

# The first derivative is a quadratic spline (order k = 3), and the formula to evaluate it
# on a knot is quite simple.
# This should give the same result as f(node.i, 0.0, Derivative(1)), but slightly faster.
function _derivative_at_node(
        ::Derivative{1}, ::CubicSplineMethod, f::ClosedFilament, node::AtNode,
    )
    (; ts, coefs, Xoffset,) = f
    (; i,) = node
    cs = coefs.cderivs[1]  # spline coefficients associated to first derivative
    t = ntuple(j -> ts[i - 2 + j], Val(3))  # = (ts[i - 1], ts[i], ts[i + 1])
    y = (
        (t[3] - t[2]) * cs[i - 1] +
        (t[2] - t[1]) * cs[i]
    ) / (t[3] - t[1])
    deperiodise_spline(y, Xoffset, ts, t[2], Val(1))  # only useful if Xoffset ≠ 0 ("infinite" / non-closed filaments)
end

# The second derivative at a node is simply equal to the corresponding spline coefficient.
function _derivative_at_node(
        ::Derivative{2}, ::CubicSplineMethod, f::ClosedFilament, node::AtNode,
    )
    f.coefs.cderivs[2][node.i]
end

function _interpolate(m::CubicSplineMethod, f::ClosedFilament, i::Int, ζ::Number, d::Derivative)
    (; ts,) = f
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]  # convert ζ ∈ [0, 1] to spline parametrisation
    _interpolate(m, f, t, d; ileft = i)
end

# Here `t_in` is in the spline parametrisation.
function _interpolate(
        ::CubicSplineMethod, f::ClosedFilament, t_in::Number, ::Derivative{n};
        ileft::Union{Nothing, Int} = nothing,
    ) where {n}
    (; ts, Xoffset,) = f
    (; cs, cderivs,) = f.coefs
    i, t = _find_knot_segment(ileft, knotlims(f), ts, t_in)
    coefs = (cs, cderivs...)[n + 1]
    ord = 4 - n
    y = evaluate_spline(coefs, ts, i, t, Val(ord))
    deperiodise_spline(y, Xoffset, ts, t_in, Val(n))  # only useful if Xoffset ≠ 0 ("infinite" / non-closed filaments)
end

# Insertion is based on the knot insertion algorithm for splines.
function _insert_node!(::CubicSplineMethod, f::ClosedFilament, i::Integer, ζ::Real)
    (; ts, Xs,) = f
    (; cs,) = f.coefs
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]  # insert knot somewhere in the middle (ζ ∈ [0, 1]; by default ζ = 0.5)
    spline_insert_knot!(cs, ts, i, t)    # note: this calls pad_periodic! on cs and ts
    Xnew = f(t; ileft = i)
    insert!(Xs, i + 1, Xnew)
    Xnew
end

function _remove_node!(::CubicSplineMethod, f::ClosedFilament, i::Integer)
    (; ts, Xs,) = f
    (; cs,) = f.coefs
    T = ts[end + 1] - ts[begin]
    popat!(cs, i)
    popat!(ts, i)
    popat!(Xs, i)
    # Avoid error if number of nodes decreased below the limit.
    # Note that `check_nodes` uses the length of f.Xs to determine this.
    if check_nodes(Bool, f)
        pad_periodic!(cs)
        pad_periodic!(ts, T)
    end
    nothing
end

# If knots were not removed but only inserted, one can pass `removed = false` to avoid
# reevaluating node positions (not needed for node insertion, since the curve is unchanged
# in that case).
function _update_after_changing_nodes!(::CubicSplineMethod, f::ClosedFilament; removed = true)
    (; Xs, ts,) = f
    (; cs,) = f.coefs
    @assert length(f) == length(Xs) == length(cs) == length(ts)  # these already have the right size
    resize!(f, length(Xs))   # resize all vectors in the filament
    check_nodes(Bool, f) || return f  # avoids error if the new number of nodes is too low
    pad_periodic!(Xs, f.Xoffset)
    # If we removed nodes, recompute spline coefficients (i.e. re-interpolate from remaining nodes).
    _update_coefficients_only!(f; only_derivatives = !removed)
    f
end
