using Test
using Random
using LinearAlgebra
using StaticArrays
using VortexFilamentEwald
using VortexFilamentEwald.Filaments

@testset "Ring" begin
    R = 3  # ring radius
    S(t) = R * SVector(1 + cospi(2t), -1 + sinpi(2t), 0)  # t ∈ [0, 1]
    Ṡ(t) = R * 2π * SVector(-sinpi(2t), cospi(2t), 0)  # derivative wrt t
    S̈(t) = R * (2π)^2 * SVector(-cospi(2t), -sinpi(2t), 0)
    S′(t) = SVector(-sinpi(2t), cospi(2t), 0)        # unit tangent (= derivative wrt the arc length ξ)
    S″(t) = SVector(-cospi(2t), -sinpi(2t), 0) ./ 3  # curvature vector (= ρ n̂ = n̂ / R)

    N = 32
    α = 1 / 100N  # amplitude of random perturbation to obtain a non-uniform node distribution

    @testset "FiniteDiffMethod(2) / HermiteInterpolation(2)" begin
        fil = @inferred Filaments.init(ClosedFilament, N, FiniteDiffMethod(2), HermiteInterpolation(2))
        ζs = collect(range(0, 1; length = N + 1)[1:N])
        rng = MersenneTwister(42)
        ζs .+= rand(rng, N) .* α
        ζs .-= ζs[begin]
        fil .= S.(ζs)
        update_coefficients!(fil)

        # Note: the actual precision depends a lot on `N` and on the perturbation amplitude `α`...
        @test isapprox(derivative(fil, 1), S′.(ζs); rtol = 1e-2)
        @test !isapprox(derivative(fil, 1), S′.(ζs); rtol = 1e-3)

        let i = 5, s = 0.3
            local t = (1 - s) * ζs[i] + s * ζs[i + 1]
            X = @inferred fil(i, s)
            Ẋ = @inferred fil(i, s, Derivative(1))
            Ẍ = @inferred fil(i, s, Derivative(2))
            X′, X″ = @inferred normalise_derivatives(Ẋ, Ẍ)
            @test isapprox(X, S(t); rtol = 1e-4)
            @test isapprox(X′, S′(t); rtol = 1e-3)
            @test isapprox(X″, S″(t); rtol = 1e-2)
        end

        # Things improve a bit when properly normalised.
        normalise_derivatives!(fil)
        @test isapprox(derivative(fil, 1), S′.(ζs); rtol = 1e-3)
        @test !isapprox(derivative(fil, 1), S′.(ζs); rtol = 1e-4)
    end

    @testset "CubicSplineMethod" begin
        fil = @inferred Filaments.init(ClosedFilament, N, CubicSplineMethod())
        ζs = collect(range(0, 1; length = N + 1)[1:N])
        rng = MersenneTwister(42)
        ζs .+= rand(rng, N) .* α
        ζs .-= ζs[begin]
        fil .= S.(ζs)
        (; ts, Xs,) = fil
        update_coefficients!(fil)

        @testset "Periodic spline knots" begin
            M = Filaments.npad(ts)
            @assert M ≥ 1
            L = ts[end + 1] - ts[begin]  # knot period (≈ total length of closed filament)
            @test isapprox(L, 2π * R; rtol = 1e-2)
            @test all(1:M) do i
                (ts[end + i] - ts[begin - 1 + i]) == (ts[end + 1 - i] - ts[begin - i]) == L
            end
        end

        @testset "Check interpolations" begin
            i = 1  # near the border, to check padding as well
            @test fil(i, 0.0) ≈ Xs[i]
            @test fil(i, 1.0) ≈ Xs[i + 1]
            @test fil(i + 1, 0.0) ≈ Xs[i + 1]
        end
    end
end
