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
    GaussLegendreQuadrature,
    quadrature

using FastGaussQuadrature: gausslegendre

"""
    AbstractQuadrature{N}

Abstract type defining an ``N``-point quadrature rule.
"""
abstract type AbstractQuadrature{N} end

function _generate(::Type{Q}) where {N, Q <: AbstractQuadrature{N}}
    f = _generator(Q)
    xs, ws = f(N)
    Xs = ntuple(i -> xs[i], Val(N))
    Ws = ntuple(i -> ws[i], Val(N))
    _change_interval(Xs, Ws)
end

# Convert quadrature rule from the [-1, 1] to the [0, 1] interval.
function _change_interval(xs, ws)
    ws′ = map(w -> w / 2, ws)
    xs′ = map(x -> (x + 1) / 2, xs)
    xs′, ws′
end

# We use a generated function to avoid runtime allocations and ensure that rules
# have static size.
"""
    quadrature(q::AbstractQuadrature{N}) -> (xs, ws)

Return ``N``-point quadrature rule valid in ``[0, 1]`` interval.

Quadrature nodes `xs` and and weights `ws` are returned as tuples.
"""
@generated quadrature(q::AbstractQuadrature) = _generate(q)

"""
    GaussLegendreQuadrature{N} <: AbstractQuadrature{N}

``N``-point Gauss–Legendre quadrature.
"""
struct GaussLegendreQuadrature{N} <: AbstractQuadrature{N} end
GaussLegendreQuadrature(N::Int) = GaussLegendreQuadrature{N}()
_generator(::Type{<:GaussLegendreQuadrature}) = gausslegendre

end
