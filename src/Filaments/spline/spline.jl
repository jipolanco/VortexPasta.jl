export
    CubicSplineMethod

using Base.Cartesian: @ntuple, @nexprs
using LinearAlgebra: ldiv!, lu!

"""
    CubicSplineMethod <: GlobalDiscretisationMethod

Represents curves using cubic splines.

In the case of closed curves, periodic cubic splines are used.
"""
struct CubicSplineMethod <: GlobalDiscretisationMethod end

# Evaluate cubic B-splines at `x` based on the knot vector `ts`.
# The `ileft` index must be such that ts[ileft] ≤ x < ts[ileft + 1].
# Note that the returned B-splines are in reversed order: (b_{i + 2}, b_{i + 1}, b_{i}, b_{i - 1}),
# where i = ileft.
function eval_cubic_bsplines(ts::AbstractVector, x::Number, ileft::Int)
    T = promote_type(eltype(ts), typeof(x))
    _evaluate_all(ts, x, Val(4), ileft, T)
end

# Evaluate cubic B-splines at ts[i].
# TODO maybe the evaluation can be further optimised knowing that we're
# evaluating right on a knot? (but it's already very fast...)
Base.@propagate_inbounds function eval_cubic_bsplines(ts::AbstractVector, i::Int)
    x = ts[i]
    bs = eval_cubic_bsplines(ts, x, i) :: NTuple{4}
    @assert iszero(first(bs))  # because we're evaluating right on the knot
    Base.tail(bs)  # return only non-zero B-splines (b_{i + 1}, b_{i}, b_{i - 1})
end

# Construct and solve cyclic tridiagonal linear system `Ax = y`.
# Uses the algorithm described in
# [Wikipedia](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants)
# based on the Sherman–Morrison formula.
function solve_cubic_spline_coefficients!(
        xs::AbstractVector, ts::PaddedVector{3}, ys::AbstractVector,
    )
    solve_cubic_spline_coefficients_slow!(xs, ts, ys)
end

# Slow version for testing only
function solve_cubic_spline_coefficients_slow!(xs, ts, ys)
    T = eltype(ts)
    n = length(ts)

    A = zeros(T, n, n)

    βs₁ = eval_cubic_bsplines(ts, 1)
    a₁ = βs₁[3]  # β₀(t₁) -- out of bands value

    βsₙ = eval_cubic_bsplines(ts, n)
    cₙ = βsₙ[1]  # βₙ₊₁(tₙ)

    ac = a₁ * cₙ
    γ = sqrt(ac)

    let bs = βs₁  # (β₂(t₁), β₁(t₁), β₀(t₁))
        A[1, 1] = bs[2] - γ  # β₁(t₁)
        A[1, 2] = bs[1]      # β₂(t₁)
    end

    for i ∈ 2:(n - 1)
        bs = eval_cubic_bsplines(ts, i)
        A[i, i - 1] = bs[3]
        A[i, i] = bs[2]
        A[i, i + 1] = bs[1]
    end

    let bs = βsₙ
        A[n, n - 1] = bs[3]
        A[n, n] = bs[2] - ac / γ
    end

    us = similar(A, n)
    fill!(us, 0)
    us[1] = γ
    us[n] = cₙ

    v₁ = 1
    vₙ = a₁ / γ

    F = lu!(A)
    ldiv!(F, us)  # us <- A \ us
    ldiv!(xs, F, ys)

    vy = v₁ * xs[1] + vₙ * xs[n]  # note: this may be a Vec{3}...
    vq = v₁ * us[1] + vₙ * us[n]  # ...while this is always a scalar

    α = vy ./ (1 + vq)
    for i ∈ eachindex(xs)
        @inbounds xs[i] = xs[i] - us[i] * α
    end

    xs
end

# TODO remove
function solve_cubic_spline_coefficients_slower!(xs, ts, ys)
    T = eltype(ts)
    n = length(ts)

    A = zeros(T, n, n)

    βs₁ = eval_cubic_bsplines(ts, 1)
    βsₙ = eval_cubic_bsplines(ts, n)

    let bs = βs₁  # (β₂(t₁), β₁(t₁), β₀(t₁))
        @assert sum(bs) ≈ 1
        A[1, 1] = bs[2]  # β₁(t₁)
        A[1, 2] = bs[1]  # β₂(t₁)
        A[1, n] = bs[3]  # β₃(t₁)
    end

    for i ∈ 2:(n - 1)
        bs = eval_cubic_bsplines(ts, i)
        @assert sum(bs) ≈ 1
        A[i, i - 1] = bs[3]
        A[i, i] = bs[2]
        A[i, i + 1] = bs[1]
    end

    let bs = βsₙ
        @assert sum(bs) ≈ 1
        A[n, n - 1] = bs[3]
        A[n, n] = bs[2]
        A[n, 1] = bs[1]
    end

    F = lu!(A)
    ldiv!(xs, F, ys)

    xs
end

# ============================================================================ #
# B-spline and spline evaluation code below adapted from BSplineKit.jl
# ============================================================================ #

@generated function _evaluate_all(
        ts::AbstractVector, x::Number, ::Val{k},
        ileft::Int, ::Type{T};
    ) where {k, T}
    @assert k ≥ 1
    ex = quote
        bs_1 = (one(T),)
    end
    for q ∈ 2:k
        bp = Symbol(:bs_, q - 1)
        bq = Symbol(:bs_, q)
        ex = quote
            $ex
            Δs = @ntuple(
                $(q - 1),
                j -> @inbounds($T(_knotdiff(x, ts, ileft - j + 1, $q - 1))),
            )
            $bq = _evaluate_step(Δs, $bp, Val($q))
        end
    end
    bk = Symbol(:bs_, k)
    quote
        $ex
        return $bk
    end
end

Base.@propagate_inbounds function _knotdiff(x::Number, ts::AbstractVector, i, n)
    j = i + n
    @boundscheck checkbounds(ts, i)
    @boundscheck checkbounds(ts, j)
    @inbounds ti = ts[i]
    @inbounds tj = ts[j]
    # @assert ti ≠ tj
    (x - ti) / (tj - ti)
end

@inline @generated function _evaluate_step(Δs, bp, ::Val{k}) where {k}
    ex = quote
        @inbounds b_1 = Δs[1] * bp[1]
    end
    for j = 2:(k - 1)
        bj = Symbol(:b_, j)
        ex = quote
            $ex
            @inbounds $bj = (1 - Δs[$j - 1]) * bp[$j - 1] + Δs[$j] * bp[$j]
        end
    end
    b_last = Symbol(:b_, k)
    quote
        $ex
        @inbounds $b_last = (1 - Δs[$k - 1]) * bp[$k - 1]
        @ntuple $k b
    end
end

# Note: there's a small change of convention compared to BSplineKit.jl, as here
# we assume that the coefficient cs[j] corresponds to the *largest* B-spline at x = ts[j].
# In other words, there is a shift in the coefficient indices wrt the
# BSplineKit implementation (see definition of `h` below).
@generated function spline_kernel(
        cs::AbstractVector, ts, i, x, ::Val{k},
    ) where {k}
    # Algorithm adapted from https://en.wikipedia.org/wiki/De_Boor's_algorithm
    h = k ÷ 2  # in BSplineKit.jl this would be just `k`
    ex = quote
        # We add zero to make sure that d_j doesn't change type later.
        # This is important when x is a ForwardDiff.Dual.
        # We also use broadcasting in case `cs` is a vector of StaticArrays.
        z = zero(x) + zero(eltype(ts))
        @nexprs $k j -> d_j = @inbounds z .+ cs[j + i - $h]
        T = typeof(d_1)
    end
    for r = 2:k, j = k:-1:r
        d_j = Symbol(:d_, j)
        d_p = Symbol(:d_, j - 1)
        jk = j - k
        jr = j - r
        ex = quote
            $ex
            @inbounds ti = ts[$jk + i]
            @inbounds tj = ts[$jr + i + 1]
            α = (x - ti) / (tj - ti)
            $d_j = ((1 - α) * $d_p + α * $d_j) :: T
        end
    end
    d_k = Symbol(:d_, k)
    quote
        $ex
        return $d_k
    end
end

# ============================================================================ #
