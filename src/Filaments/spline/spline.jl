export
    CubicSplineMethod

using Base.Cartesian: @ntuple, @nexprs
using Bumper: Bumper

"""
    CubicSplineMethod <: GlobalDiscretisationMethod

Represents curves using cubic splines.

In the case of closed curves, periodic cubic splines are used.
"""
struct CubicSplineMethod <: GlobalDiscretisationMethod end
interpolation_method(::CubicSplineMethod) = CubicSplineMethod()  # for compatibility with local methods (finite diff)

npad(::Type{<:CubicSplineMethod}) = 3  # padding needed for cubic splines

# Cubic splines are C² at the knots.
continuity(::Type{CubicSplineMethod}) = 2

struct CubicSplineCoefs{
        Method <: CubicSplineMethod,
        N,  # number of derivatives included (usually 2)
        Points <: AbstractVector,
    } <: DiscretisationCoefs{Method, N}
    method :: Method

    # B-spline coefficients associated to the curve.
    cs :: Points

    # B-spline coefficients associated to first and second derivatives.
    cderivs :: NTuple{N, Points}
end

function init_coefficients(method::CubicSplineMethod, Xs::AbstractVector, Nderiv::Val)
    cs = similar(Xs)
    cderivs = ntuple(_ -> similar(Xs), Nderiv)
    CubicSplineCoefs(method, cs, cderivs)
end

allvectors(x::CubicSplineCoefs) = (x.cs, x.cderivs...)

## ================================================================================ ##

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
        Xoffset = zero(eltype(Xs)),
    )
    periodise_coordinates!(cs, Xs, ts, Xoffset)  # useful if Xoffset ≠ 0
    n = length(ts)
    T = eltype(ts)
    # Use Bumper to allocate buffer arrays "for free" (not managed by Julia's GC).
    buf = Bumper.default_buffer()
    Bumper.@no_escape buf begin
        bufs_thomas = (
            bc = Bumper.alloc(SVector{2, T}, buf, n),
            us = Bumper.alloc(T, buf, n),
        )
        _solve_periodic_cubic_spline_coefficients_thomas!(ts, cs, bufs_thomas)
    end
    pad_periodic!(cs)
    cs
end

function periodise_coordinates!(
        Ys::AbstractVector, Xs::AbstractVector,
        ts::AbstractVector, Xoffset::Vec3,
    )
    if Xoffset === zero(Xoffset)
        copyto!(Ys, Xs)
    else
        # In this case X is not really periodic. We define Y = X - (t / T) * Xoffset
        # which is really periodic.
        # @inbounds @assert Xoffset ≈ Xs[end + 1] - Xs[begin]
        @inbounds T = ts[end + 1] - ts[begin]  # knot period
        @inbounds for i ∈ eachindex(Ys)
            # Change of variables to have a periodic spline.
            # See also `deperiodise_spline` below, which undoes this after the spline
            # has been evaluated.
            Ys[i] = Xs[i] - (ts[i] / T) * Xoffset
        end
    end
    Ys
end

# This should be used to undo periodisation of non-periodic splines.
function deperiodise_spline(y::Vec3, Xoffset::Vec3, ts::AbstractVector, t::Number, derivative::Val)
    Xoffset === zero(Xoffset) && return y
    @inbounds T = ts[end + 1] - ts[begin]  # knot period
    _deperiodise_spline(y, Xoffset, T, t, derivative)
end

_deperiodise_spline(y::Vec3, Xoffset::Vec3, T, t, ::Val{0}) = y + (t / T) * Xoffset  # zero-th derivative
_deperiodise_spline(y::Vec3, Xoffset::Vec3, T, t, ::Val{1}) = y + (1 / T) * Xoffset  # first derivative
_deperiodise_spline(y::Vec3, Xoffset::Vec3, T, t, ::Val)    = y  # higher-order derivatives

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

# Coefficient offset for periodic splines.
# This offset is such that the non-zero region of the collocation matrix is
# centred about the diagonal.
periodic_coef_offset(::Val{k}) where {k} = k ÷ 2

# Note: there's a small change of convention compared to BSplineKit.jl, as here
# we assume that the coefficient cs[j] corresponds to the *largest* B-spline at x = ts[j].
# In other words, there is a shift in the coefficient indices wrt the
# BSplineKit implementation (see definition of `h` below).
@generated function evaluate_spline(
        cs::AbstractVector, ts, i, x, ::Val{k},
    ) where {k}
    # Algorithm adapted from https://en.wikipedia.org/wiki/De_Boor's_algorithm
    # h = k ÷ 2  # in BSplineKit.jl this would be just `k`
    h = k - periodic_coef_offset(Val(k))
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

spline_derivative!(dc::PaddedVector, cs::PaddedVector, ts::PaddedVector, ord::Val) =
    spline_derivative!(copy!(dc, cs), ts, ord)

function spline_derivative!(
        cs::PaddedVector, ts::PaddedVector, ::Val{k},
    ) where {k}
    q = k - 1
    δ = periodic_coef_offset(Val(q))
    if iseven(k)
        @inbounds for i ∈ eachindex(cs)
            # We assume the denominator is non-zero, since knots should never
            # overlap in our case (they're strictly increasing).
            cs[i] = q * (cs[i + 1] - cs[i]) / (ts[i + q - δ] - ts[i - δ])
        end
    else
        @inbounds for i in Iterators.Reverse(eachindex(cs))
            cs[i] = q * (cs[i] - cs[i - 1]) / (ts[i + q - δ] - ts[i - δ])
        end
    end
    pad_periodic!(cs)
    cs
end

# Boehm's (1980) knot insertion algorithm applied to periodic cubic splines.
# The new knot `t` should be such that ts[i] ≤ t < ts[i + 1].
function spline_insert_knot!(cs::PaddedVector, ts::PaddedVector, i::Int, t::Real)
    T = ts[end + 1] - ts[begin]
    cs_new = ntuple(Val(3)) do m
        j = i - 1 + m
        @inbounds α = (t - ts[j - 2]) / (ts[j + 1] - ts[j - 2])
        @inbounds (1 - α) * cs[j - 1] + α * cs[j]
    end
    @inbounds cs[i] = cs_new[1]
    @inbounds insert!(cs, i + 1, cs_new[2])
    @inbounds cs[i + 2] = cs_new[3]
    insert!(ts, i + 1, t)
    if i + 2 > lastindex(cs)
        pad_periodic!(FromRight(), cs)   # preserve values such as cs[end + 1] (over of cs[1])
    else
        pad_periodic!(FromCentre(), cs)  # usual padding, preserves cs[1]
    end
    pad_periodic!(FromCentre(), ts, T)
    nothing
end

# ============================================================================ #
