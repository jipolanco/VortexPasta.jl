export
    CubicSplineMethod

using Base.Cartesian: @ntuple, @nexprs

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

# Obtain coefficients `cs` of periodic cubic spline from positions `Xs` and knots `ts`.
function solve_cubic_spline_coefficients!(
        cs::PaddedVector, ts::PaddedVector{3}, Xs::PaddedVector;
        buf::PaddedVector = similar(Xs),
    )
    GC.@preserve buf begin
        bufs = _thomas_buffers(ts, buf, Val(true))
        _solve_periodic_cubic_spline_coefficients_thomas!(ts, copyto!(cs, Xs), bufs)
    end
    pad_periodic!(cs)
    cs
end

function _thomas_buffers(ts::AbstractVector, buf_in::PaddedVector, ::Val{unsafe} = Val(true)) where {unsafe}
    n = length(ts)
    T = eltype(ts)
    bdata = parent(buf_in)  # this is usually a vector of Vec3{T} of length n + 2M, where M is the padding
    Base.require_one_based_indexing(bdata)
    @assert sizeof(bdata) ≥ 3 * n * sizeof(T)
    if unsafe
        # This creates two tiny allocations (96 bytes total) due to `unsafe_wrap`.
        # Tested on Julia 1.9-rc2.
        # But things are much faster than with the safe version!
        ptr = pointer(bdata)
        ptr_bc = convert(Ptr{SVector{2, T}}, ptr)
        sizeof_bc = 2n * sizeof(T)
        ptr_us = convert(Ptr{T}, ptr + sizeof_bc)
        (;
            bc = unsafe_wrap(Array, ptr_bc, n; own = false),
            us = unsafe_wrap(Array, ptr_us, n; own = false),
        )
    else
        data = reinterpret(T, bdata)
        (;
            bc = reinterpret(SVector{2, T}, view(data, 1:2n)),
            us = reinterpret(T, view(data, (2n + 1):(3n))),
        )
    end
end

# Construct and solve cyclic tridiagonal linear system `Ax = y`.
# Uses an optimised version of the algorithm described in
# [Wikipedia](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm#Variants)
# based on the Thomas algorithm and the Sherman–Morrison formula.
function _solve_periodic_cubic_spline_coefficients_thomas!(ts, xs, bufs)
    T = eltype(ts)
    n = length(ts)
    (; bc, us,) = bufs
    @assert eltype(bc) === SVector{2, T}
    @assert eltype(us) === T
    @assert length(bc) == length(us) == n

    βs₁ = eval_cubic_bsplines(ts, 1)
    βsₙ = eval_cubic_bsplines(ts, n)

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
        eval_cubic_bsplines(ts, i)  # generates i-th row of linear system
    end

    vy = v₁ * xs[1] + vₙ * xs[n]  # note: this may be a Vec{3}...
    vq = v₁ * us[1] + vₙ * us[n]  # ...while this is always a scalar

    α = vy ./ (1 + vq)
    for i ∈ eachindex(xs)
        @inbounds xs[i] = xs[i] - us[i] * α
    end

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
# ============================================================================ #
