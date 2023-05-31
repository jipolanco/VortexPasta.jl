using Test
using VortexPasta.Filaments

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

@testset "Reconnections: ∞ curve" begin
    Az = 0.001
    curve = figure_eight_curve(; a = π / 4, Az, origin = Vec3(0, 0, 0))
    (; S, tlims,) = curve
    N = 32
    ζs = range(tlims...; length = 2N + 1)[2:2:2N]
    f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())

    ts = knots(f)
    l_min = minimum(eachindex(ts)) do i
        ts[i + 1] - ts[i]
    end

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

if @isdefined(Makie)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    plot!(ax, f; refinement = 8)
    for f ∈ fs_all
        plot!(ax, f; refinement = 8)
    end
    fig
end
