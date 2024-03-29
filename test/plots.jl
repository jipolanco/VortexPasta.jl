using CairoMakie
using LinearAlgebra: norm
using Random
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using Test

function init_trefoil_filament(N::Int)
    function S(t)
        u = 2t
        π .+ (
            sinpi(u) + 2 * sinpi(2u),
            cospi(u) - 2 * cospi(2u),
            -sinpi(3u),
        )
    end
    f = @inferred Filaments.init(S, ClosedFilament, N, CubicSplineMethod())
    f
end

@testset "Makie recipes" begin
    # Just make sure that things don't error.
    Random.seed!(42)
    @testset "Filaments" begin
        f = init_trefoil_filament(16)
        @testset "Trefoil basic" begin
            plt = plot(f; axis = (type = Axis3,))  # same as `filamentplot`
            save("trefoil_basic.png", plt)
            @test isfile("trefoil_basic.png")
        end
        @testset "Trefoil" begin
            fig = Figure()
            ax = Axis3(fig[1, 1])
            plot!(  # same as `filamentplot!`
                ax, f;
                refinement = 8,
                color = :Grey,
                tangents = true, tangentcolor = :Blue,
                curvatures = true, curvaturecolor = :Red,
                vectorpos = 0.3,
            )
            save("trefoil.png", fig)
            @test isfile("trefoil.png")
        end
        @testset "Trefoil in periodic domain" begin
            periods = (1.5π, 2π, 1π)
            fig = Figure()
            ax = Axis3(fig[1, 1]; aspect = :data)
            plot!(  # same as `filamentplot!`
                ax, f;
                refinement = 8,
                color = :Grey,
                tangents = true, tangentcolor = :Blue,
                curvatures = true, curvaturecolor = :Red,
                vectorpos = 0.3,
                arrowscale = 0.5, arrowsize = (0.1, 0.1, 0.2),
                periods,
            )
            wireframe!(ax, Rect(0, 0, 0, periods...); color = :grey)
            save("trefoil_periodic.png", fig)
            @test isfile("trefoil_periodic.png")
        end
        @testset "Trefoil with velocities" begin
            velocities = randn(Vec3{Float32}, length(nodes(f)))
            fig = Figure()
            ax = Axis3(fig[1, 1]; aspect = :data)
            plot!(  # same as `filamentplot!`
                ax, f, velocities;
                refinement = 4,
                color = :Grey,
                velocitycolor = norm.(velocities),
                colormap = :cividis,
            )
            save("trefoil_velocities.png", fig)
            @test isfile("trefoil_velocities.png")
        end
    end
end
