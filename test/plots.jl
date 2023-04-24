using CairoMakie
using LinearAlgebra: norm
using Random
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using Test

function init_trefoil_filament(N::Int)
    S(t) = (
        sinpi(t) + 2 * sinpi(2t),
        cospi(t) - 2 * cospi(2t),
        -sinpi(3t),
    )
    f = Filaments.init(ClosedFilament, N, CubicSplineMethod())
    ζs = range(0, 2; length = N + 1)[1:N]
    f .= S.(ζs)
    update_coefficients!(f)
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
