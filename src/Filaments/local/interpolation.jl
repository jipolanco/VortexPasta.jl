"""
    LocalInterpolationMethod

Abstract type defining a local interpolation method for estimating curve properties
(coordinates, derivatives) in-between interpolation points using local information.
"""
abstract type LocalInterpolationMethod end

continuity(m::Union{LocalInterpolationMethod, DiscretisationMethod}) = continuity(typeof(m))

"""
    interpolate(
        method::LocalInterpolationMethod, derivative::Derivative, t::Number,
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
