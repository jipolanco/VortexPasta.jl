using Test
using Random
using LinearAlgebra
using StaticArrays
using StructArrays: StructArray, StructVector
using ForwardDiff: ForwardDiff
using VortexPasta.Quadratures: NoQuadrature, GaussLegendre
using VortexPasta.Filaments
using VortexPasta.Filaments: discretisation_method
using VortexPasta.PredefinedCurves: define_curve, Ring, TrefoilKnot
using VortexPasta.BasicTypes: VectorOfVectors
using VortexPasta.PaddedArrays: PaddedVector

function test_filament_ring(N, method)
    f = @inferred Filaments.init(ClosedFilament, N, method)

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
    if !(method isa FourierMethod)
        # Don't randomise locations for FourierMethod, as it expects roughly equal node spacings.
        γs .+= rand(rng, N) .* α
    end
    γs .-= γs[begin]

    f .= S.(γs)
    update_coefficients!(f)
    ts = knots(f)

    # Check that the interpolation actually matches the input points.
    @testset "Interpolation" begin
        err = sum(eachindex(f)) do i
            s⃗_expected = S(γs[i])
            s⃗_interp = f(i, 0.0)
            sum(abs2, s⃗_interp - s⃗_expected)
        end / N
        @test err < 1e-30
    end

    continuity = Filaments.continuity(Filaments.interpolation_method(f))

    quads = (NoQuadrature(), GaussLegendre(1), GaussLegendre(2), GaussLegendre(3), GaussLegendre(4))
    @testset "Integrate: $quad" for quad ∈ quads
        @testset "Filament length" begin
            L_expected = 2π * R  # ring perimeter
            integrand(f, i, ζ) = norm(f(i, ζ, Derivative(1)))
            L = @inferred integrate(integrand, f, quad)
            # @show method, continuity, quad, (L - L_expected) / L_expected
            rtol = if N < 10
                0.1  # just for very low resolution cases
            elseif continuity == 0 || quad === NoQuadrature()
                2e-3
            elseif method isa FiniteDiffMethod || quad === GaussLegendre(1)
                5e-3
            elseif method isa FourierMethod
                1e-9
            else
                5e-4
            end
            @test isapprox(L, L_expected; rtol)
        end
        @testset "Integration limits" begin
            # Integrate length of a single segment
            integrand(f, i, ζ) = norm(f(i, ζ, Derivative(1)))
            i = firstindex(f) + 3
            # Note: there may be small differences because integrals over subsegments are
            # computed with the same quadrature than the full integral.
            L = @inferred integrate(integrand, f, i, quad)
            L_a = @inferred integrate(integrand, f, i, quad; limits = (0, 0.3))
            L_b = @inferred integrate(integrand, f, i, quad; limits = (0.3, 1))
            # @show method, continuity, quad, (L - (L_a + L_b)) / L
            rtol = if method isa FiniteDiffMethod && quad === GaussLegendre(1)
                0.02
            elseif N < 10
                0.01
            elseif method isa FiniteDiffMethod && quad === GaussLegendre(2)
                0.005
            elseif method isa FourierMethod
                1e-9
            else
                8e-5
            end
            @test isapprox(L, L_a + L_b; rtol)
        end
    end

    @test eltype(f) === typeof(f[begin])
    @test startswith(summary(f), "$N-element $(nameof(typeof(f)))")

    @testset "Knot periodicity" begin
        M = Filaments.npad(ts)
        @assert M ≥ 1
        L = ts[end + 1] - ts[begin]  # knot period
        if !(method isa FourierMethod)
            # In general, knot period ≈ total length of closed filament, except for
            # FourierMethod.
            rtol = (N < 10) ? 0.1 : 0.01
            @test isapprox(L, 2π * R; rtol)
        end
        @test all(1:M) do i
            (ts[end + i] - ts[begin - 1 + i]) == (ts[end + 1 - i] - ts[begin - i]) == L
        end
    end

    @testset "Redistribute nodes" begin
        g = redistribute_nodes!(copy(f))
        # Check that endpoints and knot range have not been changed
        @test g[begin] ≈ f[begin]
        @test g[end + 1] ≈ f[end + 1]
        @test knotlims(g) == knotlims(f)
        # Check that the parametrisation is equidistant (can be described by a `range`).
        @test knots(g) ≈ range(ts[begin], ts[end + 1]; length = N + 1)[1:N]
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
        let i = 5
            # Third-order derivatives were not enabled.
            @test_throws ArgumentError f[i, TorsionScalar()]
            if continuity ≥ 3
                fc = @inferred Filaments.init(ClosedFilament, nodes(f), method; nderivs = Val(3))
                @test 3 == @inferred Filaments.nderivatives(fc)
                fsim = @inferred similar(fc)  # should preserve nderivs
                @test typeof(fc) === typeof(fsim)
                @test 3 == Filaments.nderivatives(fsim)
                @test @inferred(fc[i, TorsionScalar()]) == 0  # torsion is 0 for a planar ring
                if method isa FourierMethod
                    # Hermite interpolation of 3rd derivative not currently implemented.
                    @test_throws ArgumentError fc(i, 0.3, TorsionScalar())
                else
                    # Typically the case of QuinticSplineMethod
                    @test @inferred(fc(i, 0.3, TorsionScalar())) == 0
                end
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
        Filaments.recompute_parametrisation!(fc, GaussLegendre(4))  # this function may be removed in the future...
        L = -(-(knotlims(f)...))    # lower bound for filament length
        Lc = -(-(knotlims(fc)...))  # better estimation of filament length
        @assert L > 0 && Lc > 0
        @test Lc ≥ L  # generally, Lc will be greater than L
    end

    @testset "Fold periodic" begin
        Ls = (2π, 3π, 4π)
        forig = copy(f)
        Filaments.fold_periodic!(forig, Ls)
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
        @test @inferred(first(seg)) isa @inferred(eltype(seg))
    end

    if continuity ≥ 1 && N > 10
        @testset "Refinement (ring)" begin
            @testset "Coarsening" begin
                # This criterion will basically remove points such that ρℓ = ℓ/R > 2π/M.
                # Since this is a ring, in the end we should have between M/2 and M
                # segments of length ~2π/M (actually, we seem to get M/2 + 1 = 9 segments).
                M = 16
                crit = RefineBasedOnCurvature(2π, 2π / M)
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
                curv = 0.19
                crit = RefineBasedOnCurvature(curv)
                fc = copy(f)
                while true
                    ref = Filaments.refine!(fc, crit)
                    ref == (0, 0) && break
                end
                @test length(fc) > length(f)
                # Check that all ρℓ < 0.2
                ρℓ_max = 0.0
                for i ∈ eachindex(segments(fc))
                    ℓ = norm(fc[i + 1] - fc[i])
                    ρ = fc(i, 0.5, CurvatureScalar())
                    ρℓ_max = max(ρℓ_max, ρ * ℓ)
                end
                method = @inferred discretisation_method(fc)
                factor = method isa CubicSplineMethod ? 1 : 2
                @test ρℓ_max < factor * curv  # criterion is satisfied (only roughly for FiniteDiff)
            end
        end
    end

    nothing
end

# Test broadcasting a vector of filaments with a VectorOfVectors with consistent dimensions.
# This is quite convenient in timestepping.
function test_filaments_broadcasting()
    fs = [
        Filaments.init(define_curve(Ring()), ClosedFilament, 32, CubicSplineMethod()),
        Filaments.init(define_curve(TrefoilKnot()), ClosedFilament, 48, CubicSplineMethod()),
    ]
    vs = @inferred VectorOfVectors(map(similar ∘ nodes, fs))
    vv = @. vs + 3 * vs
    vf = @. vs + 3 * fs
    fv = @. fs + 3 * vs
    @test vv isa VectorOfVectors
    @test vf isa VectorOfVectors
    @test fv isa VectorOfVectors
    nothing
end

function test_init_from_vector_field(::Type{T} = Float32, method = CubicSplineMethod()) where {T}
    function taylor_green_velocity(x⃗::Vec3)
        x, y, z = x⃗
        sx, cx = sincos(x)
        sy, cy = sincos(y)
        sz = sin(z)
        Vec3(
            cx * sy * sz,
            -sx * cy * sz,
            0,
        )
    end
    function taylor_green_vorticity(x⃗::Vec3)
        x, y, z = x⃗
        sx, cx = sincos(x)
        sy, cy = sincos(y)
        sz, cz = sincos(z)
        Vec3(
            sx * cy * cz,
            cx * sy * cz,
            -2 * cx * cy * sz,
        )
    end
    x₀ = Vec3{T}(π/2 + π/3, π/2 + π/5, π/2 + 0.2)
    V = typeof(x₀)
    dτ::T = 0.11
    nsubsteps = 5
    redistribute = true

    @test @inferred(taylor_green_velocity(x₀)) isa V
    @test @inferred(Filaments.curl(taylor_green_velocity)(x₀)) isa V
    @test @inferred(taylor_green_vorticity(x₀)) isa V
    @test Filaments.curl(taylor_green_velocity)(x₀) === taylor_green_vorticity(x₀)  # identical results using AD

    # Directly from vorticity field
    f = @inferred Filaments.from_vector_field(
        ClosedFilament, taylor_green_vorticity, x₀, dτ, method;
        nsubsteps, redistribute,
    )

    # Alternative: obtain vorticity field using automatic differentiation (AD) of the velocity.
    ωf = @inferred Filaments.curl(taylor_green_velocity)
    g = @inferred Filaments.from_vector_field(ClosedFilament, ωf, x₀, dτ, method; nsubsteps, redistribute,)
    @test f ≈ g

    @test eltype(f) === V
    # Check type returned by interpolations
    @test @inferred(f(2, 0.1)) isa V
    @test @inferred(f(2, 0.1, Derivative(1))) isa V
    @test @inferred(f(2, 0.1, Derivative(2))) isa V

    err_allowed = if T === Float32
        method === QuinticSplineMethod() ? 4f-5 :
        method === CubicSplineMethod() ? 2f-4 :
        method === FiniteDiffMethod() ? 4f-4 :
        nothing
    elseif T === Float64
        method === QuinticSplineMethod() ? 3e-6 :
        method === CubicSplineMethod() ? 4e-5 :
        method === FiniteDiffMethod() ? 4e-4 :
        nothing
    end

    # Check that the filament is actually tangent to the vector field to a very good
    # approximation.
    let err = @inferred Filaments.distance_to_field(taylor_green_vorticity, f)
        # @show method, err
        @test err isa T
        @test err < err_allowed
    end

    let err = @inferred Filaments.distance_to_field(ωf, g)
        # @show method, err
        @test err isa T
        @test err < err_allowed
    end

    nothing
end

# Test nodes backed by a StructVector{Vec3{T}}.
function test_nodes_structvector(method = CubicSplineMethod())
    N = 16
    θs_pi = range(0, 2; length = N + 1)[1:N]
    T = Float32
    xs = StructArray([Vec3{T}(cospi(t), sinpi(t), 0) for t ∈ θs_pi])
    f = @inferred Filaments.init(ClosedFilament, xs, method)
    @test nodes(f) isa PaddedVector{M, Vec3{T}} where {M}
    @test parent(nodes(f)) isa StructVector{Vec3{T}}
    let g = @inferred similar(f)
        @test parent(nodes(g)) isa StructVector{Vec3{T}}
    end
    let g = @inferred similar(f, Float64)
        @test parent(nodes(g)) isa StructVector{Vec3{Float64}}
    end
    nothing
end

@testset "Filaments" begin
    @testset "Single ring" begin
        N = 32
        methods = (
            "FiniteDiff(2) / Hermite(2)" => (N, FiniteDiffMethod(2, HermiteInterpolation(2))),
            "FiniteDiff(2) / Hermite(1)" => (N, FiniteDiffMethod(2, HermiteInterpolation(1))),
            "FiniteDiff(2) / Hermite(0)" => (N, FiniteDiffMethod(2, HermiteInterpolation(0))),
            "CubicSpline" => (N, CubicSplineMethod()),
            "QuinticSpline" => (N, QuinticSplineMethod()),
            "FourierMethod" => (N, FourierMethod()),
            # These are edge cases which use a different linear solver:
            "QuinticSpline (N = 6)" => (6, QuinticSplineMethod()),
            "QuinticSpline (N = 5)" => (5, QuinticSplineMethod()),
        )
        @testset "$label" for (label, args) ∈ methods
            test_filament_ring(args...)
        end
    end
    @testset "From vector field" begin
        methods = (QuinticSplineMethod(), CubicSplineMethod(), FiniteDiffMethod())
        @testset "$method" for method ∈ methods
            test_init_from_vector_field(Float32, method)
            test_init_from_vector_field(Float64, method)
        end
    end
    @testset "StructVector nodes" begin
        methods = (QuinticSplineMethod(), CubicSplineMethod(), FiniteDiffMethod())
        @testset "$method" for method ∈ methods
            test_nodes_structvector(method)
        end
    end
    @testset "Broadcasting with VectorOfVectors" begin
        test_filaments_broadcasting()
    end
end
