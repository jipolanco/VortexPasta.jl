# Note: `only_derivatives` is not used, it's just there for compatibility with splines.
function _update_coefficients_only!(
        ::FiniteDiffMethod, f::ClosedFilament;
        only_derivatives = false,
    )
    compute_coefficients!(f.coefs, f.Xs, f.ts; Xoffset = f.Xoffset)
    f
end

function _derivative_at_node(
        ::Derivative{n}, ::FiniteDiffMethod, f::ClosedFilament, node::AtNode,
    ) where {n}
    (; cs, cderivs,) = f.coefs
    coefs = (cs, cderivs...)
    n ≤ length(cderivs) || throw(ArgumentError("FiniteDiffMethod supports up to 2 derivatives only"))
    coefs[n + 1][node.i]
end

# Calls Hermite interpolation functions
_interpolate(m::FiniteDiffMethod, args...; kws...) =
    _interpolate(interpolation_method(m), args...; kws...)
