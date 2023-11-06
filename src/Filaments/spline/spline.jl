export
    CubicSplineMethod,
    QuinticSplineMethod

using Base.Cartesian: @ntuple, @nexprs
using Bumper: Bumper, @no_escape, @alloc

struct SplineMethod{Order} <: GlobalDiscretisationMethod end

order(::Type{<:SplineMethod{Order}}) where {Order} = Order
order(m::SplineMethod) = order(typeof(m))

# For instance, cubic splines (order = 4) are C² at the knots.
continuity(::Type{M}) where {M <: SplineMethod} = order(M) - 2

npad(::Type{M}) where {M <: SplineMethod} = order(M) - 1  # e.g. npad = 3 for cubic splines

interpolation_method(m::SplineMethod) = m  # for compatibility with local methods (finite diff)

## ================================================================================ ##

"""
    CubicSplineMethod <: GlobalDiscretisationMethod

Represents curves using cubic splines.

In the case of closed curves, periodic cubic splines are used.
"""
const CubicSplineMethod = SplineMethod{4}

include("cubic.jl")

## ================================================================================ ##

"""
    QuinticSplineMethod <: GlobalDiscretisationMethod

Represents curves using quintic splines.

A quintic spline is made of polynomials of degree 5 and has global continuity ``C^4``.

In the case of closed curves, periodic quintic splines are used.
"""
const QuinticSplineMethod = SplineMethod{6}

# The linear system solver for quintic splines requires at least 7 nodes.
minimum_nodes(::QuinticSplineMethod) = 7

include("banded.jl")  # banded matrix support
include("quintic.jl")

## ================================================================================ ##

struct SplineCoefs{
        Method <: SplineMethod,
        N,  # number of derivatives included (usually 2)
        Points <: AbstractVector,
    } <: DiscretisationCoefs{Method, N}
    method :: Method

    # B-spline coefficients associated to the curve.
    cs :: Points

    # B-spline coefficients associated to first and second derivatives.
    cderivs :: NTuple{N, Points}
end

function init_coefficients(method::SplineMethod, Xs::AbstractVector, Nderiv::Val)
    cs = similar(Xs)
    cderivs = ntuple(_ -> similar(Xs), Nderiv)
    SplineCoefs(method, cs, cderivs)
end

allvectors(x::SplineCoefs) = (x.cs, x.cderivs...)

## ================================================================================ ##

# Evaluate B-splines at `x` based on the knot vector `ts`.
# The `ileft` index must be such that ts[ileft] ≤ x < ts[ileft + 1].
# Note that the returned B-splines are in reversed order.
# For example, for cubic splines: (b_{i + 2}, b_{i + 1}, b_{i}, b_{i - 1}), where i = ileft.
function eval_bsplines(order::Val, ts::AbstractVector, x::Number, ileft::Int)
    T = promote_type(eltype(ts), typeof(x))
    _evaluate_all(ts, x, order, ileft, T)
end

# TODO maybe the evaluation can be further optimised knowing that we're
# evaluating right on a knot? (but it's already very fast...)
# Evaluate B-splines at ts[i].
Base.@propagate_inbounds function eval_bsplines(order::Val{k}, ts::AbstractVector, i::Int) where {k}
    x = ts[i]
    bs = eval_bsplines(order, ts, x, i) :: NTuple{k}
    @assert iszero(first(bs))  # because we're evaluating right on the knot
    Base.tail(bs)  # return only non-zero B-splines (b_{i + 1}, b_{i}, b_{i - 1})
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
            # See also `deperiodise_curve` below, which undoes this after the spline
            # has been evaluated.
            Ys[i] = Xs[i] - (ts[i] / T) * Xoffset
        end
    end
    Ys
end

# This should be used to undo periodisation of non-periodic splines.
function deperiodise_curve(y::Vec3, Xoffset::Vec3, ts::AbstractVector, t::Number, derivative::Val)
    Xoffset === zero(Xoffset) && return y
    @inbounds T = ts[end + 1] - ts[begin]  # knot period
    _deperiodise_curve(y, Xoffset, T, t, derivative)
end

_deperiodise_curve(y::Vec3, Xoffset::Vec3, T, t, ::Val{0}) = y + (t / T) * Xoffset  # zero-th derivative
_deperiodise_curve(y::Vec3, Xoffset::Vec3, T, t, ::Val{1}) = y + (1 / T) * Xoffset  # first derivative
_deperiodise_curve(y::Vec3, Xoffset::Vec3, T, t, ::Val)    = y  # higher-order derivatives

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

# Boehm's (1980) knot insertion algorithm applied to periodic splines of order `k`.
# The new knot `t` should be such that ts[i] ≤ t < ts[i + 1].
function spline_insert_knot!(::Val{k}, cs::PaddedVector, ts::PaddedVector, i::Int, t::Real) where {k}
    T = ts[end + 1] - ts[begin]
    @assert iseven(k) "odd-order splines not currently supported"
    h = k ÷ 2
    ileft = i - h + 1  # so that ileft + 1 is the first coefficient to modify
    cs_new = ntuple(Val(k - 1)) do m
        j = ileft + m
        @inbounds α = (t - ts[j - h]) / (ts[j + h - 1] - ts[j - h])
        @inbounds (1 - α) * cs[j - 1] + α * cs[j]
    end
    # Note: the modified and new coefficients are cs[(i - h + 2):(i + h)].
    for j ∈ 1:(h - 1)
        @inbounds cs[ileft + j] = cs_new[j]
    end
    @inbounds insert!(cs, i + 1, cs_new[h])
    for j ∈ 1:(h - 1)
        @inbounds cs[i + 1 + j] = cs_new[h + j]
    end
    insert!(ts, i + 1, t)
    ifirst = i - h + 2  # index of first modified coefficient
    ilast = i + h       # index of last modified coefficient
    if ifirst < firstindex(cs)
        pad_periodic!(FromLeft(), cs)    # preserve values such as cs[begin - 1] (over cs[end])
    elseif ilast > lastindex(cs)
        pad_periodic!(FromRight(), cs)   # preserve values such as cs[end + 1] (over of cs[begin])
    else
        pad_periodic!(FromCentre(), cs)  # usual padding, preserves cs[begin]
    end
    pad_periodic!(FromCentre(), ts, T)
    nothing
end

# ============================================================================ #
