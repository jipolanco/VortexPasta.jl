using Test
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping

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

@testset "Reconnections" begin
    @testset "Static: figure-eight knot" begin
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
            @test f2 == vcat(f[begin:i], f[j + 1:end])
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
            @test f2 == vcat(f[begin:i], f[j + 1:end])
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
            Nf = length(iter.fs)

            # Run until the first self-reconnection.
            if Nf == 2
                @test 1.5 < t < 1.6
                break
            end
        end

        @test length(iter.fs) == 2
    end
end

if @isdefined(Makie)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    plot!(ax, f; refinement = 8)
    for f ∈ fs_all
        plot!(ax, f; refinement = 8)
    end
    fig
end
