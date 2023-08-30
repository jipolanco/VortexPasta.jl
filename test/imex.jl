using Test
using VortexPasta.Filaments
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot
using VortexPasta.BiotSavart
using VortexPasta.Timestepping

##

function biot_savart_parameters()
    # Grid-related parameters
    Ls = (1, 1, 1) .* 2π
    Ns = (1, 1, 1) .* 32
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    α = kmax / 5
    rcut = 4 * sqrt(2) / α

    # Physical vortex parameters
    Γ = 1.2
    a = 1e-6
    Δ = 1/4  # full core

    ParamsBiotSavart(;
        Γ, a, Δ,
        α, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = FINUFFTBackend(),
        quadrature_short = GaussLegendre(3),
        quadrature_long = GaussLegendre(5),
    )
end

##

@testset "Local/non-local splitting" begin
    ##
    params_bs = @inferred biot_savart_parameters()

    f_init = Filaments.init(
        define_curve(TrefoilKnot()), ClosedFilament, 42, CubicSplineMethod(),
    )
    fs_init = [f_init]
    tspan = (0.0, 1.0)

    prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs);

    scheme = IMEXEuler()
    iter = @inferred init(prob, scheme; dt = 0.025)

    (; fs, vs, rhs!,) = iter

    ##

    # Check that slow component (non-local) + fast component (LIA) == full velocity
    @testset "Velocity splitting" begin
        vs_full = similar(vs)
        vs_slow = similar(vs)
        vs_fast = similar(vs)

        rhs!(vs_full, fs, iter.time.t, iter; component = Val(:full))
        rhs!(vs_slow, fs, iter.time.t, iter; component = Val(:slow))
        rhs!(vs_fast, fs, iter.time.t, iter; component = Val(:fast))

        @test vs ≈ vs_full  # vs is computed when the solver is initialised
        @test vs_slow + vs_fast ≈ vs_full
    end
end
