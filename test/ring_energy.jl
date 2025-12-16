using Test
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using JET: JET

function test_vortex_ring_energy()
    # Initialise ring
    # We exaggerate the size of the ring: its diameter is almost the size of the box in the
    # periodic case.
    R = 3.0
    N = 32
    method = QuinticSplineMethod()
    S = define_curve(Ring(); scale = R)
    f = Filaments.init(S, ClosedFilament, N, method)
    fs = [f]

    # Define common Biot-Savart parameters
    params_common = (;
        Γ = 2.0,
        a = 1e-6,
        Δ = 1/4,
        quadrature = GaussLegendre(4),
    )
    tspan = (0.0, 1.0)  # this is not really important

    # Periodic case
    Ngrid = 32
    L = 2π
    Ns = (1, 1, 1) .* Ngrid
    kmax = π * Ngrid / L
    β = 3
    α = kmax / 2β
    rcut = β / α
    params_periodic = ParamsBiotSavart(; params_common..., Ls = L, rcut, Ns, α)
    prob_periodic = VortexFilamentProblem(fs, tspan, params_periodic)
    iter_periodic = init(prob_periodic, RK4(); dt = 0.1)
    if VERSION ≥ v"1.12"
        # This seems to fail on 1.11 due to Threads.Atomic.
        JET.@test_opt ignored_modules=(Base,) Diagnostics.kinetic_energy(iter_periodic)
    end
    # JET.@test_call Diagnostics.kinetic_energy(iter_periodic)  # fails when using Threads.Atomic
    E_periodic = Diagnostics.kinetic_energy(iter_periodic)
    E_periodic_quad = Diagnostics.kinetic_energy(iter_periodic; quad = GaussLegendre(4))
    @test 0 == @allocated Diagnostics.kinetic_energy(iter_periodic; nthreads = 1)
    @test 0 == @allocated Diagnostics.kinetic_energy(iter_periodic; nthreads = 1, quad = GaussLegendre(4))
    @test isapprox(E_periodic, E_periodic_quad; rtol = 1e-8)  # roughly the same result
    E_periodic_wrong = @test_logs(
        (:warn, r"should only be called when working with non-periodic domains"),
        Diagnostics.kinetic_energy_nonperiodic(iter_periodic)
    )

    # Non-periodic case
    params_nonper = ParamsBiotSavart(; params_common..., Ls = Infinity(), α = Zero())
    prob_nonper = VortexFilamentProblem(fs, tspan, params_nonper)
    iter_nonper = init(prob_nonper, RK4(); dt = 0.1)
    JET.@test_opt ignored_modules=(Base,) Diagnostics.kinetic_energy_nonperiodic(iter_nonper)
    JET.@test_call Diagnostics.kinetic_energy_nonperiodic(iter_nonper)
    E_nonper_v = Diagnostics.kinetic_energy_nonperiodic(iter_nonper)
    E_nonper_quad = Diagnostics.kinetic_energy_nonperiodic(iter_nonper; quad = GaussLegendre(4))
    @test 0 == @allocated Diagnostics.kinetic_energy_nonperiodic(iter_nonper)
    # @test 0 == @allocated Diagnostics.kinetic_energy_nonperiodic(iter_nonper; quad = GaussLegendre(4))
    @test isapprox(E_nonper_v, E_nonper_quad; rtol = 1e-8)  # roughly the same result
    E_nonper_ψ = Diagnostics.kinetic_energy_from_streamfunction(iter_nonper)

    # Compare energies obtained from both cases.
    # To be able to make comparisons, we must multiply the periodic energy (per unit mass) by
    # the domain volume, to obtain an energy per unit density.
    E_periodic_V = E_periodic * prod(BiotSavart.periods(params_periodic))

    E_expected = (params_common.Γ^2 * R/2) * (log(8R / params_common.a) - (params_common.Δ + 1))  # energy per unit density
    # @show E_periodic_wrong E_periodic_V E_nonper_v E_nonper_ψ E_expected

    # E_nonper_ψ gives the expected energy.
    # @show (E_nonper_ψ - E_expected) / E_expected
    @test isapprox(E_nonper_ψ, E_expected; rtol = 1e-3)

    # E_nonper_v overestimates the energy by a factor Γ²R/2 = Γ²/(4π) * L, where L is the
    # vortex length.
    E_overestimate = E_expected + params_common.Γ^2 * R / 2
    # @show (E_nonper_v - E_overestimate) / E_overestimate
    @test isapprox(E_nonper_v, E_overestimate; rtol = 1e-4)

    # We expect the periodic energy to be slightly smaller than the non-periodic one, for two
    # reasons:
    # 1. The non-periodic energy is the kinetic energy induced by the ring in the whole space,
    #    while the periodic one is the energy induced on a single periodic cell.
    # 2. The periodic energy includes the influence of image vortices, which tend to reduce the
    #    fluid velocity in the unit cell.
    # Still, the two values should be of the same order of magnitude.
    # @show E_periodic_V / E_nonper_ψ
    @test 0.90 < E_periodic_V / E_nonper_ψ < 0.92  # the actual difference depends on R/L and other parameters

    nothing
end

@testset "Vortex ring energy" begin
    test_vortex_ring_energy()
end
