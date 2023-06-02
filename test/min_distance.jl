# Compute the minimum distance between two filament segments.

using Test
using LinearAlgebra: norm, ⋅
using VortexPasta.Filaments

@testset "Minimum distance between filament segments" begin
    R = 0.4
    x⃗₀ = Vec3(-1.0, 0.0, 0.0)
    y⃗₀ = Vec3(+1.0, 0.0, 0.0)
    d_min_exact = abs(y⃗₀[1] - x⃗₀[1]) - 2R  # minimum distance between the circles

    N = 28
    ζs = range(0, 2; length = 3N + 1)[3:3:3N]
    xs = [x⃗₀ + R * Vec3(+cospi(t), sinpi(t), zero(t)) for t ∈ ζs]  # put the first segment on the rightmost extreme
    ys = [y⃗₀ + R * Vec3(-cospi(t), zero(t), sinpi(t)) for t ∈ ζs]  # put the first segment on the leftmost extreme
    fx = Filaments.init(ClosedFilament, xs, CubicSplineMethod())
    fy = Filaments.init(ClosedFilament, ys, CubicSplineMethod())

    # The filaments are constructed such that the following is true:
    @test isapprox(fx(0, 1/3), x⃗₀ + Vec3(R, 0, 0); rtol = 1e-5)
    @test isapprox(fy(0, 1/3), y⃗₀ - Vec3(R, 0, 0); rtol = 1e-5)

    # Case 1: the minimum distance is given by points in the interior of both segments.
    let i = 0, j = 0
        (; ζx, ζy, d⃗, x⃗, y⃗, d⃗, p⃗,) = Filaments.find_min_distance(fx, fy, i, j)
        @test d⃗ == x⃗ - y⃗ + p⃗
        @test all(iszero, p⃗)  # no periodicity -> p⃗ (periodic offset) is always zero
        @test isapprox(ζx, 1/3; rtol = 1e-3)
        @test isapprox(ζy, 1/3; rtol = 1e-3)
        @test isapprox(norm(d⃗), d_min_exact; rtol = 1e-5)
    end

    # Case 2: the optimal locations include one filament node.
    let i = 1, j = 0
        (; ζx, ζy, d⃗,) = Filaments.find_min_distance(fx, fy, i, j)
        @test ζx == 0  # minimum is at node fx[i]
        @test isapprox(ζy, 1/3; rtol = 1e-3)
        @test norm(d⃗) < norm(fx[i] - fy[j])
    end

    # Case 3: the optimal locations are both filament nodes.
    let i = 1, j = 1
        (; ζx, ζy, d⃗,) = Filaments.find_min_distance(fx, fy, i, j)
        @test ζx == 0  # minimum is at node fx[i]
        @test ζy == 0  # minimum is at node fy[j]
    end

    # This can be useful for visualising things
    if @isdefined(Makie)
        fig = Figure()
        ax = Axis3(fig[1, 1]; aspect = :data)
        filamentplot!(ax, fx; refinement = 8)
        filamentplot!(ax, fy; refinement = 8)
        for i ∈ 0:4
            data = Filaments.find_min_distance(fx, fy, i, i)
            lines!(ax, [data.x⃗, data.y⃗]; color = :grey)
        end
        fig
    end
end
