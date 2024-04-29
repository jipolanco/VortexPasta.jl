"""
    Quadratures

Module defining quadrature rules for numerical integration along filaments.

Quadrature rules are determined using the
[FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl)
package.

See the [Wikipedia page on Gaussian
quadratures](https://en.wikipedia.org/wiki/Gaussian_quadrature) for more
details.

For convenience, quadrature rules in this module are defined on the ``[0, 1]``
interval, as opposed to the standard ``[-1, 1]`` interval.
"""
module Quadratures

export
    NoQuadrature,
    GaussLegendre,
    quadrature

using FastGaussQuadrature: gausslegendre

"""
    AbstractQuadrature

Abstract type defining a quadrature rule.
"""
abstract type AbstractQuadrature end

"""
    StaticSizeQuadrature{N} <: AbstractQuadrature

Abstract type defining an ``N``-point quadrature rule.

These quadrature rules can be created basically for free.
"""
abstract type StaticSizeQuadrature{N} <: AbstractQuadrature end

"""
    PreallocatedQuadrature{T <: AbstractFloat} <: AbstractQuadrature

Abstract type defining a preallocated quadrature with element type `T`.

This is generally used for adaptive quadratures.
"""
abstract type PreallocatedQuadrature{T <: AbstractFloat} <: AbstractQuadrature end

# Convert quadrature rule from the [-1, 1] to the [0, 1] interval.
function _change_interval(xs, ws)
    ws′ = map(w -> w / 2, ws)
    xs′ = map(x -> (x + 1) / 2, xs)
    xs′, ws′
end

# This is for compatibility with PreallocatedQuadrature.
Base.convert(::Type{T}, q::StaticSizeQuadrature) where {T <: AbstractFloat} = q

"""
    Base.length(q::StaticSizeQuadrature{N}) -> N

Return the number of points associated to the quadrature.
"""
Base.length(q::StaticSizeQuadrature) = length(typeof(q))
Base.length(::Type{<:StaticSizeQuadrature{N}}) where {N} = N

"""
    GaussLegendre(N) <: StaticSizeQuadrature{N}

``N``-point Gauss–Legendre quadrature.
"""
struct GaussLegendre{N} <: StaticSizeQuadrature{N} end
GaussLegendre(N::Int) = GaussLegendre{N}()
_generator(::Type{<:GaussLegendre}) = gausslegendre

"""
    NoQuadrature() <: StaticSizeQuadrature{1}

Inexpensive 1-point quadrature rule.

When integrating, evaluates approximated values at the segment midpoint.
However, unlike a 1-point Gauss–Legendre quadrature, midpoint values are interpolated using
a basic linear interpolation.
In other words, segments are assumed to be straight.
"""
struct NoQuadrature <: StaticSizeQuadrature{1} end

"""
    Quadratures.integrate(f::Function, quad::AbstractQuadrature, lims::NTuple{2,T})
    Quadratures.integrate(f::Function, quad::AbstractQuadrature, ::Type{T})

Integrate `f(x)` using the chosen quadrature.

There are two variants:

- The first one requires passing the limits `(a, b)`, which should be floats of the desired
  accuracy;

- the second one assumes the default limits `(0, 1)`, but requires setting the float type `T`.
  It may avoid some operations when the limits are `(0, 1)`.

In both cases `T` must be a subtype of `AbstractFloat` (e.g. `Float64`).

Note that `T` is simply ignored for adaptive quadratures such as [`AdaptiveTanhSinh`](@ref),
which preallocate the quadrature coefficients with a possibly different element type `T′`
(which can be chosen when creating the quadrature).
"""
function integrate end

# Use default limits (0, 1).
function integrate(
        f::F, quad::StaticSizeQuadrature, ::Type{T},
    ) where {F <: Function, T <: AbstractFloat}
    xs, ws = quadrature(T, quad)
    vs = map(f, xs)
    u = zero(first(vs))
    for j ∈ eachindex(ws, vs)
        u = u + ws[j] * vs[j]
    end
    u
end

# Use different limits.
function integrate(
        f::F, quad::StaticSizeQuadrature, lims::NTuple{2,T},
    ) where {F <: Function, T <: AbstractFloat}
    xs, ws = quadrature(T, quad)
    a, b = lims
    δ = b - a
    ys = map(x -> a + δ * x, xs)
    vs = map(f, ys)
    u = zero(first(vs))
    for j ∈ eachindex(ws, vs)
        u = u + ws[j] * vs[j]
    end
    δ * u
end

include("tanh_sinh.jl")

## ================================================================================ ##
# The following functions must be at the end of the file to avoid world age issues
# (in particular during precompilation).

function _generate(::Type{T}, ::Type{Q}) where {N, T <: AbstractFloat, Q <: StaticSizeQuadrature{N}}
    f = _generator(Q)
    xs, ws = f(N)
    Xs = ntuple(i -> T(xs[i]), Val(N))
    Ws = ntuple(i -> T(ws[i]), Val(N))
    _change_interval(Xs, Ws)
end

# We use a generated function to avoid runtime allocations and ensure that rules
# have static size.
"""
    quadrature([T = Float64], q::StaticSizeQuadrature{N}) -> (xs, ws)

Return ``N``-point quadrature rule valid in ``[0, 1]`` interval.

Quadrature nodes `xs` and weights `ws` are returned as tuples.

This only works for fixed-size (non-adaptive) quadrature rules such as [`GaussLegendre`](@ref).
"""
@generated quadrature(::Type{T}, q::StaticSizeQuadrature) where {T <: AbstractFloat} =
    _generate(T, q)

quadrature(q::StaticSizeQuadrature) = quadrature(Float64, q)

end
