using Test
using LinearAlgebra: norm
using VortexPasta.Filaments

function trefoil_curve(; R, origin,)
    tlims = (0, 2)
    S(t) = origin .+ R .* Vec3(
        sinpi(t) + 2 * sinpi(2t),
        cospi(t) - 2 * cospi(2t),
        -sinpi(3t),
    )
    (; tlims, S, R, origin,)
end

@testset "Refinement: trefoil knot" begin
    methods = (
        FiniteDiffMethod(2),
        CubicSplineMethod(),
    )
    N = 64
    curve = trefoil_curve(; R = π / 4, origin = π .* (1, 1, 1))
    ζs = range(curve.tlims...; length = N + 1)[1:N]
    Xs_base = curve.S.(ζs)

    @testset "Method: $method" for method ∈ methods
        f = Filaments.init(ClosedFilament, copy(Xs_base), method)
        l_min_orig, l_max_orig = extrema(eachindex(segments(f))) do i
            norm(f[i + 1] - f[i])
        end

        @testset "RefineBasedOnSegmentLength" begin
            fc = copy(f)
            # Note: the max/min factor must be ≥ 1.5.
            # Otherwise we get an infinite loop:
            # 1. insert → now we're too fine
            # 2. remove → now we're too coarse
            # 3. insert → now we're too fine
            # 4. etc...
            crit = RefineBasedOnSegmentLength(2 * l_min_orig, 3 * l_min_orig)
            niter = 0
            while true  # do multiple passes if needed, until we don't need to refine anymore
                niter += 1
                ret = Filaments.refine!(fc, crit)
                all(iszero, ret) && break
                niter == 10 && break
            end
            @test niter < 10  # shouldn't need so many iterations (infinite loop?)

            l_min, l_max = extrema(eachindex(segments(fc))) do i
                # Use knots as a proxy (and approximation) for segment length.
                # This is also assumed in refine! for performance reasons.
                ts = knots(fc)
                ts[i + 1] - ts[i]
            end
            @test l_min ≥ crit.ℓ_min
            @test l_max ≤ crit.ℓ_max
        end

        @testset "RefineBasedOnCurvature" begin
            fc = copy(f)
            crit = RefineBasedOnCurvature(0.35, 0.35 / 2.2)
            niter = 0
            while true  # do multiple passes if needed, until we don't need to refine anymore
                niter += 1
                ret = Filaments.refine!(fc, crit)
                all(iszero, ret) && break
                niter == 10 && break
            end
            @test niter < 10  # shouldn't need so many iterations (infinite loop?)
            ρl_min, ρl_max = extrema(eachindex(segments(fc))) do i
                # Use knots as a proxy (and approximation) for segment length.
                # This is also assumed in refine! for performance reasons.
                ts = knots(fc)
                ℓ = ts[i + 1] - ts[i]
                ρ = (fc[i, CurvatureScalar()] + fc[i + 1, CurvatureScalar()]) / 2
                ρ * ℓ
            end
            @test ρl_min ≥ crit.ρℓ_min
            @test ρl_max ≤ crit.ρℓ_max
        end
    end
end
