# Specialisations of _derivative_at_node for cubic splines (faster than the generic
# implementation).
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
    deperiodise_curve(y, Xoffset, ts, t[2], Val(1))  # only useful if Xoffset ≠ 0 ("infinite" / non-closed filaments)
end

# The second derivative at a node is simply equal to the corresponding spline coefficient.
function _derivative_at_node(
        ::Derivative{2}, ::CubicSplineMethod, f::ClosedFilament, node::AtNode,
    )
    f.coefs.cderivs[2][node.i]
end

# Obtain coefficients `cs` of periodic spline from positions `Xs` and knots `ts`.
function solve_spline_coefficients!(
        ::Val{4}, cs::PaddedVector{M}, ts::PaddedVector{M}, Xs::PaddedVector;
        Xoffset = zero(eltype(Xs)),
    ) where {M}
    periodise_coordinates!(cs, Xs, ts, Xoffset)  # useful if Xoffset ≠ 0
    n = length(ts)
    T = eltype(ts)
    # Use Bumper to allocate buffer arrays "for free" (not managed by Julia's GC).
    buf = Bumper.default_buffer()
    @no_escape buf begin
        bufs_thomas = (
            bc = @alloc(SVector{2, T}, n),
            us = @alloc(T, n),
        )
        _solve_periodic_cubic_spline_coefficients_thomas!(cs, ts, bufs_thomas)
    end
    pad_periodic!(cs)
    cs
end

# Construct and solve cyclic tridiagonal linear system `Ax = y`.
# Uses an optimised version of the algorithm described in
# [Wikipedia](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants)
# based on the Thomas algorithm and the Sherman–Morrison formula.
function _solve_periodic_cubic_spline_coefficients_thomas!(xs, ts, bufs)
    T = eltype(ts)
    n = length(ts)
    (; bc, us,) = bufs
    @assert eltype(bc) === SVector{2, T}
    @assert eltype(us) === T
    @assert length(bc) == length(us) == n

    order = Val(4)
    βs₁ = eval_bsplines(order, ts, 1)
    βsₙ = eval_bsplines(order, ts, n)

    a₁ = βs₁[3]  # β₀(t₁)
    cₙ = βsₙ[1]  # βₙ₊₁(tₙ)
    ac = a₁ * cₙ

    # The γ coefficient is kind of arbitrary.
    # We choose the value such that the first and last elements of the main
    # diagonal (`bs`) are perturbed equally.
    γ = sqrt(ac)

    c₁ = βs₁[1]
    b₁ = βs₁[2] - γ
    bₙ = βsₙ[2] - ac / γ
    aₙ = βsₙ[3]

    fill!(us, 0)
    us[1] = γ
    us[n] = cₙ

    v₁ = 1
    vₙ = a₁ / γ

    solve_thomas!(bc, (xs, us), (b₁, c₁), (aₙ, bₙ)) do i
        @inline
        eval_bsplines(order, ts, i)  # generates i-th row of linear system
    end

    vy = v₁ * xs[1] + vₙ * xs[n]  # note: this may be a Vec{3}...
    vq = v₁ * us[1] + vₙ * us[n]  # ...while this is always a scalar

    α = vy ./ (1 + vq)
    for i ∈ eachindex(xs)
        @inbounds xs[i] = xs[i] - us[i] * α
    end

    xs
end

# Simultaneously solve M tridiagonal linear systems using Thomas algorithm.
# Optimised version interleaving construction and solution of the linear
# system, avoiding allocation of the subdiagonal `a`.
@fastmath function solve_thomas!(
        generate_row::F,                    # row generation function
        bc::AbstractVector{<:SVector{2}},   # buffers for diagonal and superdiagonal
        ds::Tuple{Vararg{AbstractVector}},  # right-hand sides, then solutions
        (b₁, c₁),                           # first row
        (aₙ, bₙ),                           # last row
    ) where {F <: Function}
    Base.require_one_based_indexing(bc)
    foreach(Base.require_one_based_indexing, ds)
    n = length(bc)
    @assert all(x -> length(x) == n, ds)
    bprev = b₁
    cprev = c₁
    bc[1] = (bprev, cprev)
    @inbounds for i ∈ 2:(n - 1)
        ci, bi, ai = generate_row(i)
        w = ai / bprev
        bprev = bi - w * cprev
        cprev = ci
        bc[i] = (bprev, cprev)
        foreach(ds) do d
            @inbounds d[i] = d[i] - w * d[i - 1]
        end
    end
    @inbounds let i = n
        bi, ai = bₙ, aₙ
        w = ai / bprev
        bc[i] = (bi - w * cprev, zero(cprev))
        foreach(ds) do d
            @inbounds d[i] = d[i] - w * d[i - 1]
        end
    end
    let bn = bc[n][1]
        foreach(ds) do d
            @inbounds d[n] = d[n] / bn
        end
    end
    for i ∈ (n - 1):-1:1
        bi, ci = bc[i]
        foreach(ds) do d
            @inbounds d[i] = (d[i] - ci * d[i + 1]) / bi
        end
    end
    ds
end

