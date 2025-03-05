using Test
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.PredefinedCurves
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using JET: JET
using KernelAbstractions: KernelAbstractions as KA  # for JET only
using LinearAlgebra: norm
using UnicodePlots: UnicodePlots, lineplot, lineplot!

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

function test_ring_stretching()
    T = Float32
    params = generate_biot_savart_parameters(T)
    # S = define_curve(TrefoilKnot(); translate = π, scale = π / 4)
    S = define_curve(Ring(); translate = π, scale = π / 4)
    N = 48
    f = Filaments.init(S, ClosedFilament, N, QuinticSplineMethod())
    fs_init = [f]
    tspan = (0.0, 1.0)
    prob = VortexFilamentProblem(fs_init, tspan, params)

    times = T[]
    energy = T[]
    vortex_length = T[]

    function callback(iter)
        local (; fs, t, dt, nstep,) = iter
        if nstep == 0
            foreach(empty!, (times, energy, vortex_length))
        end
        local quad = GaussLegendre(3)
        local E = Diagnostics.kinetic_energy(iter; quad)
        local Lvort = Diagnostics.filament_length(iter; quad)
        # if iter.stretching_velocity !== nothing
        #     write_vtkhdf("stretching_$nstep.vtkhdf", fs; refinement = 4, periods = iter.prob.p.Ls) do io
        #         io["velocity"] = iter.vs
        #     end
        # end
        push!(times, t)
        push!(energy, E)
        push!(vortex_length, Lvort)
    end

    δ = minimum(node_distance, prob.fs)
    d_crit = 0.75 * δ
    reconnect = ReconnectBasedOnDistance(d_crit)
    refinement = RefineBasedOnSegmentLength(0.75 * δ)

    dt = BiotSavart.kelvin_wave_period(params, δ) / 2
    scheme = RK4()

    # For v(ρ) = 1 / (τ * ρ), the length should exponentially increase as L(t) = L₀ * exp(t / τ).
    stretching_velocity = let L = params.Ls[1], ρ₀ = 1/L, τ = 5.0
        # ρ -> (1 - exp(-ρ / ρ₀)) / (τ * ρ)
        ρ -> -expm1(-ρ / ρ₀) / (τ * ρ)
    end
    iter = @inferred init(
        prob, scheme;
        dt, refinement, reconnect,
        stretching_velocity,
        adaptivity = NoAdaptivity(),
        callback,
    )
    @time solve!(iter)

    enable_jet = get(ENV, "JULIA_ENABLE_JET_KA_TESTS", "false") ∈ ("true", "1")  # enable JET tests involving KA kernels
    if enable_jet
        JET.@test_opt ignored_modules=(Base, KA, Base.IteratorsMD) step!(iter)
    end

    # Analytical solution
    L_expected = @. vortex_length[1] * exp(times / stretching_velocity.τ)
    # @show norm(L_expected - vortex_length) / norm(vortex_length)
    @test isapprox(vortex_length, L_expected; rtol = 1e-4)

    let plt = lineplot(
            times, energy ./ energy[1];
            title = "Ring stretching", xlabel = "Time", ylabel = "Relative change", name = "Energy",
            ylim = (1.0, 1.25)
        )
        lineplot!(plt, times, vortex_length ./ vortex_length[1]; name = "Length")
        local (; τ,) = stretching_velocity
        lineplot!(plt, times, L_expected ./ L_expected[1]; name = "Exponential")
        println(plt)
    end

    nothing
end

@testset "Ring stretching" begin
    test_ring_stretching()
end
