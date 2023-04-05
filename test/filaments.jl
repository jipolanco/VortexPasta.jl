using Test
using Random
using LinearAlgebra
using StaticArrays
using VortexFilamentEwald
using VortexFilamentEwald.Filaments

##

@testset "Ring" begin
    S(t) = 3 * SVector(1 + cospi(2t), -1 + sinpi(2t), 0)  # t ∈ [0, 1]
    Ṡ(t) = 3 * 2π * SVector(-sinpi(2t), cospi(2t), 0)  # derivative wrt t
    S̈(t) = 3 * (2π)^2 * SVector(-cospi(2t), -sinpi(2t), 0)
    S′(t) = SVector(-sinpi(2t), cospi(2t), 0)        # unit tangent (= derivative wrt the arc length ξ)
    S″(t) = SVector(-cospi(2t), -sinpi(2t), 0) ./ 3  # curvature vector (= ρ n̂ = n̂ / R)

    N = 32
    α = 1 / 100N  # amplitude of random perturbation to obtain a non-uniform node distribution

    @testset "FiniteDiff(2) / HermiteInterpolation(2)" begin
        fil = @inferred Filaments.init(ClosedFilament, N, FiniteDiff(2))
        ts = collect(range(0, 1; length = N + 1)[1:N])
        rng = MersenneTwister(42)
        ts .+= rand(rng, N) .* α
        ts .-= ts[begin]
        fil .= S.(ts)
        estimate_derivatives!(fil)

        # Note: the actual precision depends a lot on `N` and on the perturbation amplitude `α`...
        @test isapprox(derivative(fil, 1), S′.(ts); rtol = 1e-2)
        @test !isapprox(derivative(fil, 1), S′.(ts); rtol = 1e-3)

        let i = 5, s = 0.3
            local t = (1 - s) * ts[i] + s * ts[i + 1]
            X = interpolate(HermiteInterpolation(2), fil, i, s)
            Ẋ = interpolate(HermiteInterpolation(2), fil, i, s, Derivative(1))
            Ẍ = interpolate(HermiteInterpolation(2), fil, i, s, Derivative(2))
            X′, X″ = normalise_derivatives(Ẋ, Ẍ)
            @test isapprox(X, S(t); rtol = 1e-4)
            @test isapprox(X′, S′(t); rtol = 1e-3)
            @test isapprox(X″, S″(t); rtol = 1e-2)
        end

        # Things improve a bit when properly normalised.
        normalise_derivatives!(fil)
        @test isapprox(derivative(fil, 1), S′.(ts); rtol = 1e-3)
        @test !isapprox(derivative(fil, 1), S′.(ts); rtol = 1e-4)
    end
end
