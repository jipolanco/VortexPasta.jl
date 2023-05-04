using Test
using Random
using LinearAlgebra
using StaticArrays
using ForwardDiff: ForwardDiff
using VortexPasta.Filaments

function test_filament_ring(f)
    R = 3  # ring radius
    S(t) = R * SVector(1 + cospi(2t), -1 + sinpi(2t), 0)  # t ∈ [0, 1]
    Ṡ(t) = R * 2π * SVector(-sinpi(2t), cospi(2t), 0)  # derivative wrt t
    S̈(t) = R * (2π)^2 * SVector(-cospi(2t), -sinpi(2t), 0)
    S′(t) = SVector(-sinpi(2t), cospi(2t), 0)        # unit tangent (= derivative wrt the arc length ξ)
    S″(t) = SVector(-cospi(2t), -sinpi(2t), 0) ./ 3  # curvature vector (= ρ n̂ = n̂ / R)

    N = length(f)
    α = 1 / 4N  # amplitude of random perturbation to obtain a non-uniform node distribution

    γs = collect(range(0, 1; length = N + 1)[1:N])  # some arbitrary parametrisation
    rng = MersenneTwister(42)
    γs .+= rand(rng, N) .* α
    γs .-= γs[begin]

    f .= S.(γs)
    update_coefficients!(f)
    ts = knots(f)

    @test eltype(f) === typeof(f[begin])
    @test startswith(summary(f), "$N-element $(nameof(typeof(f)))")

    continuity = Filaments.continuity(Filaments.interpolation_method(f))

    @testset "Knot periodicity" begin
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
        @test f(i, 0.0) == f(ts[i]) ≈ f[i, Derivative(0)] == f[i]
        @test f(i, 1.0) == f(ts[i + 1]) ≈ f[i + 1, Derivative(0)] == f[i + 1]
        @test f(i + 1, 0.0) == f(ts[i + 1]) ≈ f[i + 1, Derivative(0)] == f[i + 1]
    end

    @testset "Derivatives at nodes" begin
        i = 1
        if continuity ≥ 1
            @test f(i, 0.0, Derivative(1)) ≈ f[i, Derivative(1)]
            @test f(i, 1.0, Derivative(1)) ≈ f[i + 1, Derivative(1)]
        end
        if continuity ≥ 2  # not the case for HermiteInterpolation(1)
            @test f(i, 0.0, Derivative(2)) ≈ f[i, Derivative(2)]
            @test f(i, 1.0, Derivative(2)) ≈ f[i + 1, Derivative(2)]
        end

        # Obtain derivatives at all nodes
        Ẋs = getindex.(Ref(f), eachindex(f), Derivative(1))
        Ẍs = getindex.(Ref(f), eachindex(f), Derivative(2))
        @test length(Ẋs) == length(Ẍs) == length(f)
        normalise_derivatives!(Ẋs, Ẍs)
        @test all(x -> norm(x) ≈ 1, Ẋs)  # unit tangents
        @test all(splat((u, v) -> 1 - u ⋅ v ≈ 1), zip(Ẋs, Ẍs))  # orthogonality with curvature vectors
    end

    @testset "Derivatives" begin
        let i = 5, ζ = 0.3
            local t = (1 - ζ) * ts[i] + ζ * ts[i + 1]
            X = @inferred f(i, ζ)
            Ẋ = @inferred f(i, ζ, Derivative(1))
            Ẍ = @inferred f(i, ζ, Derivative(2))
            @test X ≈ @inferred f(t, Derivative(0))
            @test Ẋ ≈ @inferred f(t, Derivative(1))
            @test Ẍ ≈ @inferred f(t, Derivative(2))
            @test Ẋ ≈ ForwardDiff.derivative(f, t)
            @test Ẍ ≈ ForwardDiff.derivative(
                s -> ForwardDiff.derivative(f, s),
                t,
            )
            @testset "Normalisation" begin
                X′, X″ = @inferred normalise_derivatives(Ẋ, Ẍ)
                @test norm(X′) ≈ 1
                @test 1 - X′ ⋅ X″ ≈ 1   # orthogonality
                if continuity ≥ 1  # don't run test with HermiteInterpolation{0} (→ curvature = 0)
                    @test isapprox(norm(X″), 1 / R; rtol = 0.1)  # ring curvature (some methods are more accurate than others...)
                end
                @test X′ ≈ f(i, ζ, UnitTangent())
                @test X″ ≈ f(i, ζ, CurvatureVector())
            end
        end
    end

    @testset "Geometric quantities" begin
        let i = 5, ζ = 0.3
            t̂ = @inferred f(i, ζ, UnitTangent())
            ρ⃗ = @inferred f(i, ζ, CurvatureVector())
            ρ = @inferred f(i, ζ, CurvatureScalar())
            b⃗ = @inferred f(i, ζ, CurvatureBinormal())
            @test norm(t̂) ≈ 1
            @test norm(ρ⃗) ≈ ρ
            @test norm(b⃗) ≈ ρ
            @test abs(t̂ ⋅ ρ⃗) < 1e-12
            @test t̂ × ρ⃗ ≈ b⃗
            let t = (1 - ζ) * ts[i] + ζ * ts[i + 1]
                @test ρ⃗ ≈ f(t, CurvatureVector())
            end
        end
        let i = 5
            if continuity ≥ 2
                @test f[i, CurvatureVector()] ≈ f(i, 0.0, CurvatureVector())
                @test f[i + 1, CurvatureVector()] ≈ f(i, 1.0, CurvatureVector())
            end
        end
    end

    @testset "copy" begin
        fc = @inferred copy(f)
        @test typeof(f) === typeof(fc)
        Fil = typeof(f)
        # Check that all numerical data has been copied
        @test all(1:fieldcount(Fil)) do i
            getfield(f, i) == getfield(fc, i)
        end
    end

    @testset "Recompute parametrisation" begin
        fc = copy(f)
        Filaments.recompute_parametrisation!(fc)  # this function may be removed in the future...
        L = -(-(knotlims(f)...))    # lower bound for filament length
        Lc = -(-(knotlims(fc)...))  # better estimation of filament length
        @assert L > 0 && Lc > 0
        @test Lc ≥ L  # generally, Lc will be greater than L
    end

    @testset "Fold periodic" begin
        Ls = (2π, 3π, 4π)
        forig = Filaments.fold_periodic!(copy(f), Ls)
        fc = copy(forig)
        for i ∈ eachindex(fc)
            fc[i] += 2 * Vec3(Ls)
        end
        @test !(nodes(fc) ≈ nodes(forig))
        Filaments.fold_periodic!(fc, Ls)
        @test nodes(fc) ≈ nodes(forig)
    end
end

@testset "Ring" begin
    N = 32
    filaments = (
        "FiniteDiff(2) / Hermite(2)" => @inferred(Filaments.init(ClosedFilament, N, FiniteDiffMethod(2), HermiteInterpolation(2))),
        "FiniteDiff(2) / Hermite(1)" => @inferred(Filaments.init(ClosedFilament, N, FiniteDiffMethod(2), HermiteInterpolation(1))),
        "FiniteDiff(2) / Hermite(0)" => @inferred(Filaments.init(ClosedFilament, N, FiniteDiffMethod(2), HermiteInterpolation(0))),
        "CubicSpline" => @inferred(Filaments.init(ClosedFilament, N, CubicSplineMethod())),
    )
    @testset "$s" for (s, f) ∈ filaments
        test_filament_ring(f)
    end
end
