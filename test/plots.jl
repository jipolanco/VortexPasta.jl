using CairoMakie
using VortexFilamentEwald.Filaments
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
    @testset "Filaments" begin
        f = init_trefoil_filament(16)
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
        let
            plt = plot(f; axis = (type = Axis3,))  # same as `filamentplot`
            save("trefoil_basic.png", plt)
            @test isfile("trefoil_basic.png")
        end
    end
end
