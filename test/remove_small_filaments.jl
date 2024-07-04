# Timestepping: check that filaments with less than N nodes are removed at the end of an
# iteration. Here N depends on the actual discretisation method (e.g. N = 5 for
# QuinticSplineMethod).

using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.PredefinedCurves: Ring, define_curve
using VortexPasta.Timestepping
using Test

function generate_biot_savart_parameters(::Type{T}) where {T}
    Γ = 1.0
    a = 1e-6
    Δ = 1/4
    L = 2π
    Ls = (L, L, L)
    β = 3.0
    rcut = L / 2
    α = β / rcut
    kmax = 2α * β
    M = ceil(Int, kmax * L / π)
    Ns = (1, 1, 1) .* M
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short = NaiveShortRangeBackend(),
        backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4)),
        quadrature = GaussLegendre(3),
    )
end

function test_filament_removal()
    T = Float32
    params = generate_biot_savart_parameters(T)
    method = QuinticSplineMethod()

    curves = (
        define_curve(Ring(); translate = (0, 0, T(π)/4)),
        define_curve(Ring(); translate = (0, 0, T(π)/2)),
    )
    Nps = (Filaments.minimum_nodes(method), 10)

    fs = @inferred map(curves, Nps) do S, Np
        Filaments.init(S, ClosedFilament{T}, Np, method)
    end |> collect

    # Remove nodes from first filament until we're below the required number of nodes (5 for
    # QuinticSplineMethod)
    let f = fs[1]
        @assert length(f) == Filaments.minimum_nodes(method)
        @test Filaments.check_nodes(Bool, f) === true  # filament has enough nodes
        Filaments.remove_node!(f, 1)
        Filaments.update_after_changing_nodes!(f)        # doesn't fail
        @test Filaments.check_nodes(Bool, f) === false   # not enough nodes
    end

    prob = @inferred VortexFilamentProblem(fs, (0, 1), params)
    iter = @inferred init(prob, RK4(); dt = 0.1)

    # Check that the "invalid" filament was removed from the filament list.
    @test length(iter.fs) == length(iter.vs) == length(iter.ψs) == 1
    @test iter.fs[1] == prob.fs[2]
    @test length(iter.vs[1]) == length(iter.ψs[1]) == length(iter.fs[1])
    @test iter.stats.filaments_removed_count == 1
    @test iter.stats.filaments_removed_length > 0

    nothing
end

@testset "Timestepping: filament removal" begin
    test_filament_removal()
end
