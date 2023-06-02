using Test
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using LinearAlgebra: norm

# Create a curve which resembles an 8 (or ∞).
# This specific curve is called the lemniscate of Bernouilli (https://en.wikipedia.org/wiki/Lemniscate_of_Bernoulli)
# We allow a perturbation in the 3rd direction (Az) so that the curve doesn't exactly cross itself.
function figure_eight_curve(; a::Real = 1, origin::Vec3 = π * Vec3(1, 1, 1), Az = 0)
    c = a / sqrt(2)  # focal distance (`a` is the half width)
    # Note: the crossings are located at t = -0.5 and 0.5.
    tlims = (-1, 1)
    S(t) = let
        s, c = sincospi(t)
        factor = a / (1 + s^2)
        z = Az * s
        origin + Vec3(factor * c, factor * s * c, z)
    end
    (; tlims, S, a, c, Az, origin,)
end

function trefoil_knot_curve(; a::Real, origin::Vec3 = π * Vec3(1, 1, 1))
    tlims = (-1, 1)
    S(t) = origin + a * Vec3(
        sinpi(t) + 2 * sinpi(2t),
        cospi(t) - 2 * cospi(2t),
        -sinpi(3t),
    )
    (; tlims, S, a, origin,)
end

function ring_curve(; R::Real, origin::Vec3 = π * Vec3(1, 1, 1), sign = +1)
    tlims = (0, 2)
    S(t) = origin + R * Vec3(cospi(sign * t), sinpi(sign * t), zero(t))
    (; tlims, S, origin, R,)
end

@testset "Reconnections" begin
    @testset "Static self-reconnection: figure-eight knot" begin
        Az = 0.001
        curve = figure_eight_curve(; a = π / 4, Az, origin = Vec3(0, 0, 0))
        (; S, tlims,) = curve
        N = 32
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())

        l_min = minimum_knot_increment(f)

        @testset "Filaments.split!" begin
            i = length(f) ÷ 4
            j = 3i
            f1, f2 = @inferred Filaments.split!(copy(f), i, j)
            update_coefficients!.((f1, f2))
            @test length(f1) == length(f2) == length(f) ÷ 2
            @test f1 == f[i + 1:j]
            @test f2 == vcat(f[j + 1:end], f[begin:i])
        end

        @testset "reconnect_self!" begin
            crit = @inferred BasedOnDistance(l_min / 2)
            fc = copy(f)
            fs_all = [fc]
            Filaments.reconnect_self!(crit, fc, fs_all)
            @test length(fs_all) == 2  # filament split into two!
            f1, f2 = fs_all[1], fs_all[2]

            # Check that the reconnection happened at the right place
            i = length(f) ÷ 4
            j = 3i
            @test length(f1) == length(f2) == length(f) ÷ 2
            @test f1 == f[i + 1:j]
            @test f2 == vcat(f[j + 1:end], f[begin:i])
        end
    end

    @testset "Dynamic: trefoil knot" begin
        curve = trefoil_knot_curve(; a = π / 4,)
        (; S, tlims,) = curve
        N = 64
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]

        f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
        fs = [f]
        l_min = minimum_knot_increment(fs)

        params_bs = let
            Ls = (1, 1, 1) .* 2π
            Ns = (1, 1, 1) .* 32
            kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
            α = kmax / 5
            rcut = 4 / α
            ParamsBiotSavart(;
                Γ = 2.0, α, a = 1e-6, Δ = 1/4, rcut, Ls, Ns,
                backend_short = NaiveShortRangeBackend(),
                backend_long = FINUFFTBackend(),
                quadrature_short = GaussLegendreQuadrature(4),
                quadrature_long = GaussLegendreQuadrature(4),
            )
        end

        tspan = (0.0, 1.0)  # ignored
        prob = @inferred VortexFilamentProblem(fs, tspan, params_bs)
        iter = @inferred init(
            prob, RK4();
            dt = 1.0,  # will be changed by the adaptivity
            dtmin = 1e-4,
            refinement = BasedOnCurvature(π / 8; ℓ_min = l_min / 2),
            reconnect = BasedOnDistance(l_min),
            adaptivity = BasedOnSegmentLength(1.0),
        )

        # TODO check that energy decreases?
        for n = 1:200
            Timestepping.step!(iter)
            (; t,) = iter.time
            Nf = length(fs)

            # Run until the first self-reconnection.
            if Nf == 2
                @test 1.5 < t < 1.6
                break
            end
        end

        @test length(fs) == 2
    end

    @testset "Static reconnection: antiparallel rings" begin
        # Test two nearly overlapping rings with the same sign (they should reconnect,
        # since they're antiparallel at the reconnection point).
        R = π / 4
        ϵ = 0.02R
        rings = (
            ring_curve(; R, origin = Vec3(π - (R + ϵ), π, π), sign = +1),
            ring_curve(; R, origin = Vec3(π + (R + ϵ), π, π), sign = +1),
        )
        d_min = norm(rings[2].S(1.0) - rings[1].S(0.0))  # minimum distance between the rings
        @assert d_min ≈ 2ϵ

        fs_orig = map(rings) do ring
            (; S, tlims,) = ring
            N = 32
            ζs = range(tlims...; length = 2N + 1)[2:2:2N]
            Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
        end
        l_min = minimum(minimum_knot_increment, fs_orig)
        d_min_nodes = let (f, g) = fs_orig
            minimum(Iterators.product(eachindex(f), eachindex(g))) do (i, j)
                norm(f[i] - g[j])
            end
        end
        crit = BasedOnDistance(1.2 * d_min_nodes)

        @testset "Filaments.reconnect_other!" begin
            f, g = copy.(fs_orig)
            @test Filaments.reconnect_other!(crit, f, g) == true  # reconnection happened!
            @test length(f) == sum(length, fs_orig)
            @test maximum_knot_increment(f) < 2 * l_min  # check that there are no crazy jumps
        end

        @testset "Filaments.reconnect!" begin
            fs = collect(copy.(fs_orig))
            Filaments.reconnect!(crit, fs)
            @test length(fs) == 1
            @test length(fs[1]) == sum(length, fs_orig)
            @test maximum_knot_increment(fs) < 2 * l_min  # check that there are no crazy jumps
        end
    end

    @testset "Static reconnection: nearly overlapping rings" begin
        # Test two nearly overlapping rings (with some small offset in z) which should reconnect at two points.
        # We should end up with two separate vortices.
        # Internally, this happens in two steps:
        #   (1) the two vortices merge at one of the reconnection points, and
        #   (2) the resulting vortex splits into two at the other reconnection point.
        R = π / 4
        ϵ = 0.02R
        rings = (
            ring_curve(; R, origin = Vec3(π - R / sqrt(2), π, π + ϵ), sign = +1),
            ring_curve(; R, origin = Vec3(π + R / sqrt(2), π, π - ϵ), sign = +1),
        )

        fs_orig = map(rings) do ring
            (; S, tlims,) = ring
            N = 32
            ζs = range(tlims...; length = 2N + 1)[2:2:2N]
            Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
        end
        d_min_nodes = let (f, g) = fs_orig
            minimum(Iterators.product(eachindex(f), eachindex(g))) do (i, j)
                norm(f[i] - g[j])
            end
        end
        crit = BasedOnDistance(1.2 * d_min_nodes)

        fs = collect(copy.(fs_orig))
        Filaments.reconnect!(crit, fs)
        @test length(fs) == 2
        @test fs[1] != fs_orig[1]  # reconnection happened
        @test fs[2] != fs_orig[2]  # reconnection happened
        @test length(fs[1]) ≠ length(fs[2])  # the resulting vortices have different sizes
    end
end

# This can be useful for playing around with the tests.
if @isdefined(Makie)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    for f ∈ fs_orig
        plot!(ax, f; refinement = 8, linestyle = :dash)
    end
    for f ∈ fs
        plot!(ax, f; refinement = 8, linestyle = :solid)
    end
    fig
end
