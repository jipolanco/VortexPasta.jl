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
    coefs[n + 1][node.i]
end

_interpolate(m::FiniteDiffMethod, args...; kws...) =
    _interpolate(interpolation_method(m), args...; kws...)

function _interpolate(
        m::HermiteInterpolation, f::ClosedFilament, t_in::Number, d::Derivative = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    )
    (; ts,) = f
    i, t = _find_knot_segment(ileft, knotlims(f), ts, t_in)
    ζ = (t - ts[i]) / (ts[i + 1] - ts[i])
    y = _interpolate(m, f, i, ζ, d)
    _deperiodise_finitediff(d, y, f, t, t_in)
end

function _interpolate(
        method::HermiteInterpolation{M}, f::ClosedFilament,
        i::Int, t::Number, deriv::Derivative{N},
    ) where {M, N}
    (; cs, cderivs,) = f.coefs
    ts = knots(f)
    checkbounds(f, i)
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
    α * interpolate(method, deriv, t, values_i...)
end

function _deperiodise_finitediff(::Derivative{0}, y, f::ClosedFilament, t, t_in)
    (; Xoffset,) = f
    Xoffset === zero(Xoffset) && return y
    t == t_in && return y
    ta, tb = knotlims(f)
    T = tb - ta
    dt = t_in - t  # expected to be a multiple of T
    y + dt / T * Xoffset
end

_deperiodise_finitediff(::Derivative, y, args...) = y  # derivatives (n ≥ 1): shift not needed

function _insert_node!(::FiniteDiffMethod, f::ClosedFilament, i::Integer, ζ::Real)
    (; Xs,) = f
    (; cs, cderivs,) = f.coefs
    @assert length(cderivs) == 2

    s⃗ = f(i, ζ)  # new point to be added
    insert!(Xs, i + 1, s⃗)  # note: this uses derivatives at nodes (Hermite interpolations), so make sure they are up to date!
    insert!(cs, i + 1, s⃗)

    # We also insert new derivatives in case we need to perform Hermite interpolations later
    # (for instance if we insert other nodes).
    # Derivatives will be recomputed later when calling `update_after_changing_nodes!`.
    s⃗′ = f(i, ζ, Derivative(1))
    s⃗″ = f(i, ζ, Derivative(2))
    insert!(cderivs[1], i + 1, s⃗′)
    insert!(cderivs[2], i + 1, s⃗″)

    s⃗
end

function _remove_node!(::FiniteDiffMethod, f::ClosedFilament, i::Integer)
    (; ts, Xs,) = f
    # No need to modify interpolation coefficients, since they will be updated later in
    # `update_after_changing_nodes!`. This assumes that we won't insert nodes after removing
    # nodes (i.e. insertions are done before removals).
    popat!(ts, i)
    popat!(Xs, i)
end

# The `removed` argument is just there for compatibility with splines.
function _update_after_changing_nodes!(::FiniteDiffMethod, f::ClosedFilament; removed = true)
    (; Xs,) = f
    resize!(f, length(Xs))   # resize all vectors in the filament
    if check_nodes(Bool, f)  # avoids error if the new number of nodes is too low
        pad_periodic!(Xs, f.Xoffset)
        update_coefficients!(f)
    end
    f
end
