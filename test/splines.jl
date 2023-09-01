# Test evaluation of some curve properties exactly on spline knots when using cubic splines.
# The idea is to be able to evaluate the LIA term (proportional to s⃗′ × s⃗″) directly from
# the cubic spline coefficients `cs`.

using Test
using VortexPasta.Filaments
using LinearAlgebra: ×
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot

@testset "Spline evaluation on knots" begin
    f = @inferred Filaments.init(
        define_curve(TrefoilKnot()), ClosedFilament, 42, CubicSplineMethod(),
    )

    # Evaluate quantities on a chosen knot.
    i = 4
    s⃗′ = @inferred f[i, Derivative(1)]
    s⃗″ = @inferred f[i, Derivative(2)]

    # Compute quantities directly from spline coefficients.
    (; ts, cs, cderivs,) = f
    cs′, cs″ = cderivs

    # A linear spline (order k = 2) evaluated on a knot is simply equal to the spline
    # coefficient associated to that knot.
    @test s⃗″ == cs″[i]

    # Evaluation of a quadratic spline (order k = 3) on a knot.
    @test s⃗′ == (
        (ts[i + 1] - ts[i]) * cs′[i - 1] +
        (ts[i] - ts[i - 1]) * cs′[i]
    ) / (ts[i + 1] - ts[i - 1])

    # Evaluation of s⃗′ × s⃗″ from spline coefficients of s⃗′ (which is a quadratic spline).
    @test s⃗′ × s⃗″ ≈ 2 * (cs′[i - 1] × cs′[i]) / (ts[i + 1] - ts[i - 1])

    # Check relation between coefficients of s⃗ (cubic spline) and of s⃗′ (quadratic spline).
    @test cs′[i] ≈ 3 * (cs[i + 1] - cs[i]) / (ts[i + 2] - ts[i - 1])

    # Evaluation of s⃗′ × s⃗″ directly from spline coefficients of s⃗ (cubic spline).
    @test s⃗′ × s⃗″ ≈ 18 * (
        (cs[i] - cs[i - 1]) × (cs[i + 1] - cs[i])
    ) / (
        (ts[i + 1] - ts[i - 1]) *
        (ts[i + 1] - ts[i - 2]) *
        (ts[i + 2] - ts[i - 1])
    )
end
