function _update_coefficients_only!(
        method::SplineMethod, f::ClosedFilament;
        only_derivatives = false,
    )
    (; coefs, Xs, ts, Xoffset,) = f
    @assert method === coefs.method
    compute_coefficients!(coefs, Xs, ts; Xoffset, only_derivatives)
    f
end

# Compute coefficients of the derivatives of a spline.
# Note that the derivative of a spline of order `k` is a spline of order `k - 1`.
@inline function _compute_derivative_coefs!(::Val{k}, cderivs::NTuple, ck, ts) where {k}
    # Compute derivative coefficients `cd` from spline coefficients `ck`.
    cd = first(cderivs)
    spline_derivative!(cd, ck, ts, Val(k))
    _compute_derivative_coefs!(Val(k - 1), Base.tail(cderivs), cd, ts)  # compute next derivative from `cd`
end

@inline _compute_derivative_coefs!(::Val, ::Tuple{}, args...) = nothing

# Generic implementation, simply evaluates spline on the node.
function _derivative_at_node(
        d::Derivative, ::SplineMethod, f::ClosedFilament, node::AtNode,
    )
    f(node.i, 0.0, d)
end

function _interpolate(m::SplineMethod, f::ClosedFilament, i::Int, ζ::Number, d::Derivative)
    (; ts,) = f
    t = @inbounds (1 - ζ) * ts[i] + ζ * ts[i + 1]  # convert ζ ∈ [0, 1] to spline parametrisation
    _interpolate(m, f, t, d; ileft = i)
end

# Here `t_in` is in the spline parametrisation.
function _interpolate(
        method::SplineMethod, f::ClosedFilament, t_in::Number, deriv::Derivative{n};
        ileft::Union{Nothing, Int} = nothing,
    ) where {n}
    (; ts, Xoffset,) = f
    y = evaluate(method, f.coefs, f.ts, t_in, deriv; ileft)
    deperiodise_curve(y, Xoffset, ts, t_in, Val(n))  # only useful if Xoffset ≠ 0 ("infinite" / non-closed filaments)
end

# Insertion is based on the knot insertion algorithm for splines.
function _insert_node!(method::SplineMethod, f::ClosedFilament, i::Integer, ζ::Real)
    (; ts, Xs,) = f
    (; cs,) = f.coefs
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]  # insert knot somewhere in the middle (ζ ∈ [0, 1]; by default ζ = 0.5)
    k = order(method)
    spline_insert_knot!(Val(k), cs, ts, i, t)    # note: this calls pad_periodic! on cs and ts
    Xnew = f(t; ileft = i)
    insert!(Xs, i + 1, Xnew)
    Xnew
end

function _remove_node!(::SplineMethod, f::ClosedFilament, i::Integer)
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
function _update_after_changing_nodes!(::SplineMethod, f::ClosedFilament; removed = true)
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
