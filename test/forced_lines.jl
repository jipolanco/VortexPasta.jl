# Test initially straight infinite lines subject to velocity forcing.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm
using Random: Random
using VortexPasta.PredefinedCurves: define_curve, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using UnicodePlots: lineplot, lineplot!
using JET: JET
using FINUFFT: FINUFFT

function filaments_from_functions(funcs::NTuple{N, F}, args...) where {N, F <: Function}
    isempty(funcs) && return ()
    S = first(funcs)
    f = Filaments.init(S, args...)
    f_next = filaments_from_functions(Base.tail(funcs), args...)
    (f, f_next...)
end

function test_forced_lines(
        ::Type{T} = Float64;
        method = QuinticSplineMethod(),
        scheme = RK4(),
    ) where {T}
    ## Initial condition
    L = 2π
    p = PeriodicLine()  # unperturbed straight line
    funcs = let L = L
        (
            define_curve(p; scale = L, translate = (1L/4, 1L/4, L/2), orientation = +1),
            define_curve(p; scale = L, translate = (3L/4, 1L/4, L/2), orientation = -1),
            define_curve(p; scale = L, translate = (1L/4, 3L/4, L/2), orientation = -1),
            define_curve(p; scale = L, translate = (3L/4, 3L/4, L/2), orientation = +1),
        )
    end
    N = 32
    filaments = collect(filaments_from_functions(funcs, ClosedFilament{T}, N, method))

    ## Biot–Savart parameters
    Γ = 2.0
    a = 1e-6
    Δ = 1/4
    Ngrid = 16
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    Ls = (L, L, L)
    α = kmax / 5
    rcut = 5 / α

    params = ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
        # backend_short = CellListsBackend(2),
        backend_short = NaiveShortRangeBackend(),  # needed when rcut > L/2 (which is the case here, since Ngrid is small)
        backend_long = FINUFFTBackend(tol = 1e-6),
        quadrature = GaussLegendre(2),
    )

    # Simulation timestep (also appearing in the amplitude of the stochastic forcing)
    dt = BiotSavart.kelvin_wave_period(params, L/N)

    ## Forcing
    # We apply a stochastic velocity forcing.
    # The forcing changes randomly at each timestep (not the best possible implementation of
    # stochastic forcing!).
    # This is similar to the forcing used by Baggaley & Laurie (PRB 2014).
    rng = Random.Xoshiro()
    forcing_amplitude = Ref(sqrt(Γ / dt / 100))  # this has the units of a velocity [L¹T⁻¹]
    forcing_ks = T[-3, -2, -1, 1, 2, 3]  # forcing wavenumbers (in the z direction)
    forcing_ϕs = similar(forcing_ks)     # forcing phases
    forcing_as = similar(forcing_ks)     # normal-distributed amplitudes (with unit variance)

    function forcing_velocity(x⃗::Vec3, t)
        A = forcing_amplitude[]
        z = x⃗.z
        w = zero(Complex{T})
        for (k, ϕ, a) ∈ zip(forcing_ks, forcing_ϕs, forcing_as)
            w += a * cis(k * z + ϕ)
        end
        A * Vec3(real(w), imag(w), 0)
    end

    # Note: this satisfies v = ∇ × ψ
    function forcing_streamfunction(x⃗::Vec3, t)
        A = forcing_amplitude[]
        z = x⃗.z
        w = zero(Complex{T})
        for (k, ϕ, a) ∈ zip(forcing_ks, forcing_ϕs, forcing_as)
            w -= (a / k) * cis(k * z + ϕ)
        end
        A * Vec3(real(w), imag(w), 0)
    end

    ## Callback function to be called after each timestep
    times = T[]
    energy = T[]
    line_length = T[]

    function callback(iter)
        (; nstep, t,) = iter
        Tend = iter.prob.tspan[2]
        if nstep == 0  # make sure vectors are empty at the beginning of the simulation
            Random.seed!(rng, 42)  # restart the random number generator
            empty!(times)
            empty!(energy)
            empty!(line_length)
        end
        # Set random amplitudes and phases of forcing for next iteration.
        let
            Random.randn!(rng, forcing_as)  # normal distribution
            Random.rand!(rng, forcing_ϕs)
            forcing_ϕs .*= 2π  # random in [0, 2π]
        end
        E = Diagnostics.kinetic_energy_from_streamfunction(iter)
        Lvort = Diagnostics.filament_length(iter)
        push!(times, t)
        push!(energy, E)
        push!(line_length, Lvort)
        # @show nstep, t/Tend, E, Lvort/4
    end

    ## Initialise solver
    tspan = (0.0, 10 * dt)  # run a few timesteps
    prob = VortexFilamentProblem(filaments, tspan, params)

    # We first test *without* the external streamfunction, which means that energy estimates
    # will only include the energy associated to the Biot–Savart velocity (which is probably
    # the most relevant energy).
    fill!(forcing_as, 0)  # this disables forcing of the initial velocity
    fill!(forcing_ϕs, 0)
    iter = init(
        prob, scheme;
        dt,
        callback,
        external_velocity = forcing_velocity,
        # external_streamfunction = forcing_streamfunction,
    )
    @time solve!(iter)
    let
        plt = lineplot(times, line_length ./ first(line_length); title = "Forced lines", name = "Line length", xlabel = "Time", ylabel = "L/L₀")
        lineplot!(plt, times, energy ./ first(energy); name = "Energy")
        display(plt)
    end
    L_growth = last(line_length) / first(line_length) - 1
    E_growth = last(energy) / first(energy) - 1
    @show L_growth E_growth
    @test 0.004 < L_growth < 0.006
    @test 0.004 < E_growth < 0.006

    # Now test including the external_streamfunction parameter (energy values will be much
    # higher). Note that this doesn't change the dynamics, and the vortex trajectories
    # should be exactly the same as before.
    energy_prev = copy(energy)
    line_length_prev = copy(line_length)
    fill!(forcing_as, 0)
    fill!(forcing_ϕs, 0)
    iter = init(
        prob, scheme;
        dt,
        callback,
        external_velocity = forcing_velocity,
        external_streamfunction = forcing_streamfunction,
    )
    @time solve!(iter)
    @test line_length ≈ line_length_prev  # we get the exact same results
    # Energy estimates are much larger when including energy associated to the
    # forcing velocity.
    # @show energy ./ energy_prev
    @test all(10 .< energy ./ energy_prev .< 100)

    ## Inference tests
    # Check that everything is inferred (except for issues coming from FINUFFT.jl).
    # (This is only for tests and can be ignored!)
    JET.@test_opt ignored_modules=(FINUFFT, Base, Base.PCRE) init(
        prob, scheme;
        dt,
        callback,
        external_velocity = forcing_velocity,
        external_streamfunction = forcing_streamfunction,
    )
    JET.@test_opt ignored_modules=(FINUFFT, Base, Base.PCRE) step!(iter)

    nothing
end

@testset "Velocity forcing" begin
    test_forced_lines()
end
