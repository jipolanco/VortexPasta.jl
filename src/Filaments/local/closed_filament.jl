# Note: `only_derivatives` is not used, it's just there for compatibility with splines.
function _update_coefficients_only!(
        ::FiniteDiffMethod, f::ClosedFilament;
        only_derivatives = false,
    )
    (; ts, Xs,) = f
    (; cs, cderivs,) = f.coefs
    M = npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts
    method = discretisation_method(f)
    copy!(cs, Xs)  # this also copies padding (uses copyto! implementation in PaddedArrays)
    @inbounds for i ∈ eachindex(Xs)
        ℓs_i = ntuple(j -> ts[i - M + j] - ts[i - M - 1 + j], Val(2M))  # = ℓs[(i - M):(i + M - 1)]
        Xs_i = ntuple(j -> Xs[i - M - 1 + j], Val(2M + 1))  # = Xs[(i - M):(i + M)]
        coefs_dot = coefs_first_derivative(method, ℓs_i)
        coefs_ddot = coefs_second_derivative(method, ℓs_i)
        cderivs[1][i] = sum(splat(*), zip(coefs_dot, Xs_i))  # = ∑ⱼ c[j] * x⃗[j]
        cderivs[2][i] = sum(splat(*), zip(coefs_ddot, Xs_i))
    end
    # These paddings are needed for Hermite interpolations and stuff like that.
    # (In principle we just need M = 1 for two-point Hermite interpolations.)
    map(pad_periodic!, cderivs)
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
