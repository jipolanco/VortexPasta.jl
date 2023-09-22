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
        FiniteDiffMethod(),
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
            # Note: the max/min factor must be ≳ 2.
            # Otherwise we get an infinite loop:
            # 1. insert → now we're too fine
            # 2. remove → now we're too coarse
            # 3. insert → now we're too fine
            # 4. etc...
            crit = RefineBasedOnSegmentLength(2 * l_min_orig, 4 * l_min_orig)
            niter = 0

            while true  # do multiple passes if needed, until we don't need to refine anymore
                niter += 1
                ret = Filaments.refine!(fc, crit)
                all(iszero, ret) && break
                niter == 10 && break
            end
            @test niter < 10  # shouldn't need so many iterations (infinite loop?)

            l_min, l_max = extrema(segments(fc)) do seg
                norm(seg(1.0) - seg(0.0))
            end
            @test l_min ≥ crit.ℓ_min
            @test l_max ≤ crit.ℓ_max
        end

        @testset "RefineBasedOnCurvature" begin
            fc = copy(f)
            crit = RefineBasedOnCurvature(0.35, 0.35 / 2.5)
            niter = 0
            while true  # do multiple passes if needed, until we don't need to refine anymore
                niter += 1
                ret = Filaments.refine!(fc, crit)
                all(iszero, ret) && break
                niter == 10 && break
            end
            @test niter < 10  # shouldn't need so many iterations (infinite loop?)
            ρl_min, ρl_max = extrema(segments(fc)) do seg
                ℓ = norm(seg(1.0) - seg(0.0))
                ρ = (seg(0.0, CurvatureScalar()) + seg(1.0, CurvatureScalar())) / 2
                ρ * ℓ
            end
            @test ρl_min ≥ crit.ρℓ_min
            @test ρl_max ≤ crit.ρℓ_max
        end

        @testset "Knot removal" begin
            # Check that removing a node doesn't change the curve properties (continuity,
            # periodicity...). We test all nodes just in case, but nodes near the beginning
            # and end are the ones which used to fail for splines (fixed since)...
            for i ∈ eachindex(f)
                @testset "Node $i/$N" begin
                    fc = copy(f)
                    Filaments.remove_node!(fc, i)
                    Filaments.update_after_changing_nodes!(fc; removed = true)
                    @test all(i -> fc[i] ≈ fc(i, 0.0), eachindex(fc))
                    @test fc[end + 1] == fc[begin]
                    @test fc[end] == fc[begin - 1]
                    ts = knots(fc)
                    @test fc(ts[begin]) ≈ fc[begin]
                    @test fc(lastindex(fc), 1.0) ≈ fc[end + 1]
                end
            end
        end
    end
end

if @isdefined(Makie)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    # wireframe!(ax, Rect(0, 0, 0, Ls...); color = :grey, linewidth = 0.5)
    plot!(ax, f; refinement = 8, linestyle = :dash, color = (:grey, 0.5))
    let f = fc
        p = plot!(ax, f; refinement = 8)
        scatter!(ax, nodes(f).data; color = p.color)
        inds = (firstindex(f) - 2):(lastindex(f) + 2)
        for i ∈ inds
            align = i ∈ eachindex(f) ? (:left, :bottom) : (:right, :top)
            color = i ∈ eachindex(f) ? (p.color, 1.0) : (p.color, 0.6)
            text!(ax, f[i]; text = string(i), fontsize = 16, color, align)
        end
    end
    fig
end
