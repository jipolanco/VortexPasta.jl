using Test
using Random
using LinearAlgebra
using StaticArrays
using ForwardDiff: ForwardDiff
using VortexPasta.Filaments

function test_filament_ring(args)
    f = @inferred Filaments.init(ClosedFilament, args...)

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

    @testset "Segments" begin
        seg = @inferred segments(f)
        inds = @inferred eachindex(seg)
        @test inds == eachindex(nodes(f))  # true for closed filaments (same number of segments and nodes)
        @test firstindex(seg) == first(inds)
        @test lastindex(seg) == last(inds)
    end

    if continuity ≥ 1
        @testset "Refinement (ring)" begin
            @testset "Coarsening" begin
                # This criterion will basically remove points such that ρℓ = ℓ/R > 2π/M.
                # Since this is a ring, in the end we should have between M/2 and M
                # segments of length ~2π/M (actually, we seem to get M/2 + 1 = 9 segments).
                M = 16
                crit = BasedOnCurvature(2π, 2π / M)
                fc = copy(f)
                # Several iterations are needed to remove required nodes, since we never
                # remove adjacent nodes in a single pass.
                while true
                    n_add, n_rem = Filaments.refine!(fc, crit)
                    @test n_add == 0
                    @test n_rem ≥ 0
                    n_rem == 0 && break  # we're done removing nodes
                end
                @test length(fc) < length(f)
                @test M ÷ 2 ≤ length(fc) ≤ M
                # Check that all ρℓ > 2π / 2M
                ρℓ_min = Inf
                for i ∈ eachindex(segments(fc))
                    ℓ = norm(f[i + 1] - f[i])
                    ρ = f(i, 0.5, CurvatureScalar())
                    ρℓ_min = min(ρℓ_min, ρ * ℓ)
                end
                @test ρℓ_min > 2π / 4M  # criterion is (roughly) satisfied
            end
            @testset "Refining" begin
                crit = BasedOnCurvature(0.2)
                fc = copy(f)
                while true
                    Filaments.refine!(fc, crit) == (0, 0) && break
                end
                @test length(fc) > length(f)
                # Check that all ρℓ < 0.2
                ρℓ_max = 0.0
                for i ∈ eachindex(segments(fc))
                    ℓ = norm(fc[i + 1] - fc[i])
                    ρ = fc(i, 0.5, CurvatureScalar())
                    ρℓ_max = max(ρℓ_max, ρ * ℓ)
                end
                factor = fc isa ClosedSplineFilament ? 1 : 2
                @test ρℓ_max < factor * 0.2  # criterion is satisfied (only roughly for FiniteDiff)
            end
        end
    end

    nothing
end

@testset "Ring" begin
    N = 32
    methods = (
        "FiniteDiff(2) / Hermite(2)" => (N, FiniteDiffMethod(2), HermiteInterpolation(2)),
        "FiniteDiff(2) / Hermite(1)" => (N, FiniteDiffMethod(2), HermiteInterpolation(1)),
        "FiniteDiff(2) / Hermite(0)" => (N, FiniteDiffMethod(2), HermiteInterpolation(0)),
        "CubicSpline" => (N, CubicSplineMethod()),
    )
    @testset "$label" for (label, args) ∈ methods
        test_filament_ring(args)
    end
end
