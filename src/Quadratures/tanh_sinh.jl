export AdaptiveTanhSinh

"""
    AdaptiveTanhSinh([T = Float64]; nlevels = 10, rtol = sqrt(eps(T))) <: PreallocatedQuadrature{T}

Adaptive [tanh-sinh quadrature](https://en.wikipedia.org/wiki/Tanh-sinh_quadrature).

Behaves well when there are singularities at (or near) the endpoints.

It can be easily made adaptive because it can be written as a simple trapezoidal rule after
a change of variables.

# Optional arguments

- `T = Float64`: quadrature precision;

- `nlevels = 10`: maximum number of adaptivity levels. Must be `≥ 2`;

- `rtol = sqrt(eps(T))`: relative tolerance.

Computations are stopped either when the maximum adaptivity level is reached, or when the
difference between two levels falls below the relative tolerance `rtol`.

Note that the maximum number of function evaluations is `2^nlevels`, so it can make sense to
use a small number of levels (or a large tolerance) when function evaluations are expensive.
"""
struct AdaptiveTanhSinh{T <: AbstractFloat} <: PreallocatedQuadrature{T}
    xs :: Vector{T}  # locations in [0, 1] (sorted according to access order, for performance reasons)
    ws :: Vector{T}  # weights (also sorted). Note that the weights are scaled by the step Δt at the smallest level.
    w₀ :: T          # weight associated to location t = x = 0
    hmax :: T        # integration step `h` at the first level (then it's divided by 2 at each level)
    rtol :: T
    nlevels :: Int

    function AdaptiveTanhSinh(::Type{T}; nlevels::Int = 10, rtol = sqrt(eps(T))) where {T <: AbstractFloat}
        nlevels < 2 && throw(ArgumentError("`nlevels` must be ≥ 2"))
        N = 1 << nlevels  # this gives 0 if nlevels is too large (such that 2^nlevels can't be represented by an Int)
        N > 0 || throw(ArgumentError(lazy"number of levels is too large: $nlevels"))

        tmax = tanh_sinh_determine_tmax(T)
        ts = range(0, tmax; length = N + 1)[(begin + 1):end]  # we store t > 0 only (we use symmetry for t < 0)
        hmax = last(ts)  # same as tmax...

        # Permutation so that so that nodes and weights are sorted in access order.
        ps = sortperm_for_adaptive_quadrature(nlevels)

        # Compute nodes and weigths, sorted in access order.
        xs = similar(ts)
        ws = similar(ts)
        for i ∈ eachindex(ps, xs, ws)
            j = ps[i]
            t = ts[j]
            x, w = tanh_sinh_node_weight(t)
            xs[i] = x
            ws[i] = w
        end

        # Make sure we're not evaluating at the endpoints (in case of singularities...).
        # This could fail if we choose a tmax which is too large for the wanted precision T,
        # such that T(0.99999...) == T(1).
        # Also note that xs[begin] contains the largest value in xs (due to the way nodes
        # are sorted).
        @assert xs[begin] < 1
        @assert xs[begin] === prevfloat(one(T))

        _, w₀ = tanh_sinh_node_weight(zero(T))

        # Convert to [0, 1] range.
        # TODO do we need this?
        # for i ∈ eachindex(xs, ws)
        #     ws[i] = ws[i] / 2
        #     xs[i] = (xs[i] + 1) / 2
        # end

        new{T}(xs, ws, w₀, hmax, rtol, nlevels)
    end
end

AdaptiveTanhSinh(; kws...) = AdaptiveTanhSinh(Float64; kws...)

# We recompute everything just in case we're going to higher precision (e.g. Float32 -> Float64)
function Base.convert(::Type{T}, q::AdaptiveTanhSinh{S}) where {T <: AbstractFloat, S}
    (; nlevels,) = q
    if eps(T) > eps(S)  # we're going to lower precision
        rtol = max(T(q.rtol), sqrt(eps(T)))  # just in case the original tolerance is too large for T
    else
        rtol = T(q.rtol)
    end
    AdaptiveTanhSinh(T; nlevels, rtol,)
end

Base.convert(::Type{T}, q::AdaptiveTanhSinh{T}) where {T <: AbstractFloat} = q  # nothing to do

function Base.show(io::IO, q::AdaptiveTanhSinh{T}) where {T}
    (; nlevels, rtol,) = q
    print(io, lazy"AdaptiveTanhSinh($T; nlevels = $nlevels, rtol = $rtol)")
end

# Determine truncation of integral in the infinite domain.
# Since the weights decay very fast (exponentially?), it makes sense to apply a cut-off
# somewhere.
# Value found by hand such that xmax == prevfloat(1.0).
# That is, the rightmost node `xmax` is smaller than 1 by the tiniest possible amount
# allowed by the type T.
tanh_sinh_determine_tmax(::Type{Float64}) = 3.1909e0
tanh_sinh_determine_tmax(::Type{Float32}) = 2.40f0

# For t ∈ [-∞, ∞], this returns a location in [-1, 1] and its corresponding weight.
function tanh_sinh_node_weight(t)
    pi_half = oftype(t, π) / 2
    c_sinh_t = pi_half * sinh(t)
    x = tanh(c_sinh_t)
    w = (pi_half * cosh(t)) / cosh(c_sinh_t)^2
    x, w
end

function sortperm_for_adaptive_quadrature(nlevels)
    N = 1 << nlevels  # this gives 0 if nlevels is too large (such that 2^nlevels can't be represented by an Int)
    @assert N > 0
    ks = Vector{Int}(undef, N)
    l = firstindex(ks) - 1

    # Level 0
    n = N
    ks[l += 1] = n

    # Levels ≥ 1
    n = n >> 1
    m = 1
    while n ≠ 0
        for j ∈ 1:m
            k = (2j - 1) * n
            ks[l += 1] = k
        end
        n = n >> 1
        m = m << 1
    end

    ks
end

# Assume the nodes and weights are re-sorted in the order they are accessed by adaptive
# quadrature.
function integrate(f::F, quad::AdaptiveTanhSinh{T}, lims::NTuple{2,Real}) where {F <: Function, T <: AbstractFloat}
    (; xs, ws, w₀, hmax, rtol,) = quad
    Base.require_one_based_indexing(xs)
    Base.require_one_based_indexing(ws)

    N = length(xs)
    # @assert ispow2(N)
    onehalf = 1 / T(2)

    a, b = lims
    Δ = onehalf * T(b - a)
    c = onehalf * T(a + b)

    # Initial estimate: consider only central node.
    h = hmax * Δ
    Iprev = 2 * h * w₀ * f(c)

    # Second estimate
    n = unsigned(N)
    l = firstindex(xs)
    I = let
        w, x = @inbounds ws[l], xs[l]
        @inline onehalf * Iprev + h * w * (
            f(c + Δ * x) +    # rightmost node
            f(c - Δ * x)      # leftmost node
        )
    end

    # Start iterating
    n = n >>> 1
    h = onehalf * h  # divide integration step by 2 for next level
    m = unsigned(1)
    @fastmath while n ≠ 0
        Iprev = I
        I = zero(I)
        # Start by adding odd terms of the current level.
        # Note that weights must be multiplied by h (done afterwards).
        for _ ∈ 1:m
            l += 1
            w, x = @inbounds ws[l], xs[l]
            y = Δ * x
            @inline I += w * (f(c + y) + f(c - y))
        end
        I = h * I + onehalf * Iprev  # add all even terms of the current level
        Idiff_sq = sum(abs2, I - Iprev)  # should also work with SVectors
        I_rtol_sq = sum(abs2, I * rtol)
        if Idiff_sq < I_rtol_sq
            break
        end
        n = n >>> 1
        h = onehalf * h
        m = m << 1
    end

    I
end

# Default [0, 1] limits.
# Ignore `S` in this case, it's just there for compatibility with `integrate` with
# GaussLegendre or other quadratures.
integrate(f::F, quad::AdaptiveTanhSinh{T}, ::Type{S}) where {F, T, S <: AbstractFloat} =
    integrate(f, quad, (zero(T), one(T)))
