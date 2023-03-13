export
    interpolate

"""
    InterpolationMethod

Abstract type defining an interpolation method for estimating curve properties
(coordinates, derivatives) in-between interpolation points.
"""
abstract type InterpolationMethod end

include("hermite.jl")

"""
    interpolate(
        method::InterpolationMethod, [derivative = Derivative(0)], t::Number,
        values_at_interpolation_points...,
    )

Interpolate coordinates ``\\bm{X}(t)`` or their derivatives at a given
location ``t``.

Under this form, the location ``t`` must be normalised such that the interval of
interest is in ``t ∈ [0, 1]``. Note that input and output derivatives must also
be normalised accordingly.

The `values_at_interpolation_points` depend on the interpolation method:

- for [`HermiteInterpolation`](@ref), one should pass the coordinates
  ``\\bm{X}`` and the first ``M`` derivatives at the two endpoints of the interval
  (``t = 0`` and ``1``).
"""
function interpolate end

interpolate(method::InterpolationMethod, t::Number, args...) =
    interpolate(method, Derivative(0), t, args...)
