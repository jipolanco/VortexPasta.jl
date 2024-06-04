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

# This is in particular to check that knot insertion in spline implementations gives the
# right interpolation coefficients.
function check_coefficients(f)
    fc = copy(f)
    # Check that coefficients don't change after calling update_coefficients!.
    # If they do, it's because the original coefficients are wrong (they don't correctly
    # interpolate the filament nodes).
    update_coefficients!(fc; knots = knots(fc))
    @test fc.coefs.cs ≈ f.coefs.cs
    nothing
end

@testset "Refinement: trefoil knot" begin
    methods = (
        FiniteDiffMethod(),
        CubicSplineMethod(),
        QuinticSplineMethod(),
        FourierMethod(),
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

        @testset "Insert at beginning" begin
            fc = copy(f)
            i = firstindex(fc)  # insert node between nodes 1 and 2
            Filaments.insert_node!(fc, i, 0.5)
            Filaments.update_after_changing_nodes!(fc; removed = false)  # we didn't remove any nodes
            check_coefficients(fc)
        end

        @testset "Insert at end" begin
            fc = copy(f)
            i = lastindex(fc)  # insert node after last node
            Filaments.insert_node!(fc, i, 0.5)
            Filaments.update_after_changing_nodes!(fc; removed = false)  # we didn't remove any nodes
            check_coefficients(fc)
        end

        @testset "RefineBasedOnSegmentLength" begin
            # Note: the max/min factor must be ≳ 2.
            # Otherwise we get an infinite loop:
            # 1. insert → now we're too fine
            # 2. remove → now we're too coarse
            # 3. insert → now we're too fine
            # 4. etc...
            α_low = method isa FourierMethod ? 1.9 : 2.0
            crit = RefineBasedOnSegmentLength(α_low * l_min_orig, 4 * l_min_orig)

            fc = copy(f)
            niter = 0
            while true  # do multiple passes if needed, until we don't need to refine anymore
                niter += 1
                N_prev = length(fc)
                ret = Filaments.refine!(fc, crit)
                check_coefficients(fc)
                all(iszero, ret) && break
                n_add, n_rem = ret
                @test length(fc) == N_prev + n_add - n_rem
                niter == 10 && break
            end
            @test niter < 10  # shouldn't need so many iterations (infinite loop?)
            fc

            l_min, l_max = extrema(segments(fc)) do seg
                norm(seg(1.0) - seg(0.0))
            end
            @test l_min ≥ crit.ℓ_min
            @test l_max ≤ crit.ℓ_max
        end

        @testset "RefineBasedOnCurvature" begin
            fc = copy(f)
            crit = if method isa FourierMethod
                RefineBasedOnCurvature(0.2, 0.2 / 3.5)
            else
                RefineBasedOnCurvature(0.2, 0.2 / 2.5)
            end
            niter = 0
            while true  # do multiple passes if needed, until we don't need to refine anymore
                niter += 1
                ret = Filaments.refine!(fc, crit)
                check_coefficients(fc)
                all(iszero, ret) && break
                n_add, n_rem = ret
                # Make sure we insert nodes at the first iteration (since the
                # RefineBasedOnSegmentLength test mainly removes nodes).
                @assert n_add > 0 || niter > 1
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
            # Note: this is not really relevant for FourierMethod...
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
        inds = (firstindex(f) - 1):(lastindex(f) + 1)
        for i ∈ inds
            align = i ∈ eachindex(f) ? (:left, :bottom) : (:right, :top)
            color = i ∈ eachindex(f) ? (p.color, 1.0) : (p.color, 0.6)
            text!(ax, f[i]; text = string(i), fontsize = 16, color, align)
        end
    end
    fig
end
