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
    rtol :: T
    nlevels :: Int

    function AdaptiveTanhSinh(::Type{T}; nlevels::Int = 10, rtol = sqrt(eps(T))) where {T <: AbstractFloat}
        nlevels < 2 && throw(ArgumentError("`nlevels` must be ≥ 2"))
        N = 1 << nlevels  # this gives 0 if nlevels is too large (such that 2^nlevels can't be represented by an Int)
        N > 0 || throw(ArgumentError(lazy"number of levels is too large: $nlevels"))

        tmax = tanh_sinh_determine_tmax(T)
        ts = range(0, tmax; length = N + 1)[(begin + 1):end]  # we store t > 0 only (we use symmetry for t < 0)
        h = step(ts)

        # Permutation so that so that nodes and weights are sorted in access order.
        ps = sortperm_for_adaptive_quadrature(nlevels)

        # Compute nodes and weigths, sorted in access order.
        xs = similar(ts)
        ws = similar(ts)
        for i ∈ eachindex(ps, xs, ws)
            j = ps[i]
            t = ts[j]
            xs[i] = tanh_sinh_node(t)
            ws[i] = h * tanh_sinh_weight(t)  # note: weighted for step h associated to maximum level
        end

        w₀ = h * tanh_sinh_weight(zero(T))

        # Convert to [0, 1] range.
        # TODO do we need this?
        # for i ∈ eachindex(xs, ws)
        #     ws[i] = ws[i] / 2
        #     xs[i] = (xs[i] + 1) / 2
        # end

        new{T}(xs, ws, w₀, rtol, nlevels)
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
# This value is chosen since tanh_sinh_weight(3.3) ~ 1e-17.
# TODO do something that depends on the precision, based on decay of tanh_sinh_weight(t)?
# For example, use a smaller value for Float32?
tanh_sinh_determine_tmax(::Type{T}) where {T <: AbstractFloat} = 3.3

# For t ∈ [-∞, ∞], this returns a location in [-1, 1].
function tanh_sinh_node(t)
    pi_half = oftype(t, π) / 2
    tanh(pi_half * sinh(t))
end

function tanh_sinh_weight(t)
    pi_half = oftype(t, π) / 2
    (pi_half * cosh(t)) / cosh(pi_half * sinh(t))^2
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
    (; xs, ws, w₀, rtol,) = quad
    Base.require_one_based_indexing(xs)
    Base.require_one_based_indexing(ws)

    N = length(xs)
    # @assert ispow2(N)
    onehalf = 1 / T(2)

    a, b = lims
    Δ = onehalf * T(b - a)
    c = onehalf * T(a + b)

    # Initial estimate: consider only central node.
    Iprev = 2 * N * w₀ * f(c)

    # Second estimate
    n = unsigned(N)
    l = firstindex(xs)
    I = let
        w, x = @inbounds ws[l], xs[l]
        @inline onehalf * Iprev + n * w * (
            f(c + Δ * x) +    # rightmost node
            f(c - Δ * x)      # leftmost node
        )
    end

    # Start iterating
    n = n >>> 1
    m = unsigned(1)
    @fastmath while n ≠ 0
        Iprev = I
        I = zero(I)
        # Start by adding odd terms of the current level.
        # Note that weights must be multiplied by n (done afterwards).
        for _ ∈ 1:m
            l += 1
            w, x = @inbounds ws[l], xs[l]
            y = Δ * x
            @inline I += w * (f(c + y) + f(c - y))
        end
        I = n * I + onehalf * Iprev  # add all even terms of the current level
        Idiff_sq = sum(abs2, I - Iprev)  # should also work with SVectors
        I_rtol_sq = sum(abs2, I * rtol)
        if Idiff_sq < I_rtol_sq
            break
        end
        n = n >>> 1
        m = m << 1
    end

    # @assert l == N + 1

    Δ * I
end

# Default [0, 1] limits.
# Ignore `S` in this case, it's just there for compatibility with `integrate` with
# GaussLegendre or other quadratures.
integrate(f::F, quad::AdaptiveTanhSinh{T}, ::Type{S}) where {F, T, S <: AbstractFloat} =
    integrate(f, quad, (zero(T), one(T)))
