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

# ================================================================================ #

abstract type AbstractStochasticForcing{T <: AbstractFloat} end

function init_forcing!(rng::Random.AbstractRNG, f::AbstractStochasticForcing)
    Random.randn!(rng, f.ws)  # normal distribution
    f.ws .*= f.amplitude
    f
end

function generate_forcing_function(f::AbstractStochasticForcing{T}) where {T}
    function (x⃗::Vec3, t)
        z = x⃗.z
        u = zero(Complex{T})
        for (k, w) ∈ zip(f.ks, f.ws)
            u += w * cis(k * z)
        end
        Vec3(real(u), imag(u), 0)
    end
end

function generate_streamfunction_function(f::AbstractStochasticForcing{T}) where {T}
    # Note: this satisfies v = ∇ × ψ
    function (x⃗::Vec3, t)
        z = x⃗.z
        u = zero(Complex{T})
        for (k, w) ∈ zip(f.ks, f.ws)
            u += (-w / k) * cis(k * z)
        end
        Vec3(real(u), imag(u), 0)
    end
end

# ================================================================================ #

# Delta-correlated forcing (like Brownian motion)
struct BrownianForcing{T} <: AbstractStochasticForcing{T}
    dt :: T          # simulation timestep
    σ² :: T          # fixed forcing amplitude -- has units of a diffusion coefficient [L²T⁻¹]
    amplitude :: T   # = √{σ² / dt} -- has units of a velocity [LT⁻¹]
    ks :: Vector{T}  # forced wavenumbers
    ws :: Vector{Complex{T}}  # random forcing coefficients (~ Normal(0, std = amplitude))
end

correlation_time(f::BrownianForcing) = f.dt

function BrownianForcing(dt, σ², ks)
    ws = similar(ks, complex(eltype(ks)))
    fill!(ws, 0)
    amplitude = sqrt(σ² / dt)
    BrownianForcing(dt, σ², amplitude, ks, ws)
end

function update_forcing!(rng::Random.AbstractRNG, f::BrownianForcing)
    init_forcing!(rng, f)  # we just reset the coefficients to new random values
end

# ================================================================================ #

# Stochastic forcing with correlation time τ
struct OrnsteinUhlenbeckForcing{T} <: AbstractStochasticForcing{T}
    dt :: T          # simulation timestep
    σ² :: T          # fixed forcing amplitude -- has units of a diffusion coefficient [L²T⁻¹]
    amplitude :: T   # = √{σ² / dt} -- has units of a velocity [LT⁻¹]
    τ  :: T          # correlation time
    ks :: Vector{T}  # forced wavenumbers
    ws :: Vector{Complex{T}}  # random forcing coefficients
end

correlation_time(f::OrnsteinUhlenbeckForcing) = f.τ

function OrnsteinUhlenbeckForcing(dt, σ², τ_, ks)
    τ = max(τ_, dt)  # numerically, the correlation time can't be smaller than dt
    amplitude = sqrt(σ² / dt)
    ws = similar(ks, complex(eltype(ks)))
    fill!(ws, 0)
    OrnsteinUhlenbeckForcing(dt, σ², amplitude, τ, ks, ws)
end

function update_forcing!(rng::Random.AbstractRNG, f::OrnsteinUhlenbeckForcing)
    (; ws, dt, τ, amplitude,) = f
    for (i, w) ∈ pairs(ws)
        # Note: setting τ = dt gives Brownian motion
        ws[i] = (1 - dt / τ) * w + amplitude * Random.randn(rng, typeof(w)) # * dt
    end
    f
end

# ================================================================================ #

function filaments_from_functions(funcs::NTuple{N, F}, args...) where {N, F <: Function}
    isempty(funcs) && return ()
    S = first(funcs)
    f = Filaments.init(S, args...)
    f_next = filaments_from_functions(Base.tail(funcs), args...)
    (f, f_next...)
end

function generate_straight_lines(::Type{T}, L, N, method) where {T}
    p = PeriodicLine()  # unperturbed straight line
    funcs = let L = L
        (
            define_curve(p; scale = L, translate = (1L/4, 1L/4, L/2), orientation = +1),
            define_curve(p; scale = L, translate = (3L/4, 1L/4, L/2), orientation = -1),
            define_curve(p; scale = L, translate = (1L/4, 3L/4, L/2), orientation = -1),
            define_curve(p; scale = L, translate = (3L/4, 3L/4, L/2), orientation = +1),
        )
    end
    collect(filaments_from_functions(funcs, ClosedFilament{T}, N, method))
end

function generate_biot_savart_parameters()
    L = 2π
    Γ = 9.97e-4
    a = 1e-8
    Δ = 1/4
    Ngrid = 16
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    Ls = (L, L, L)
    α = kmax / 5
    rcut = 5 / α
    ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
        # backend_short = CellListsBackend(2),
        backend_short = NaiveShortRangeBackend(),  # needed when rcut > L/2 (which is the case here, since Ngrid is small)
        backend_long = FINUFFTBackend(tol = 1e-6),
        quadrature = GaussLegendre(2),
    )
end

function test_forced_lines(
        ::Type{T} = Float64;
        N = 32,
        method = QuinticSplineMethod(),
        scheme = RK4(),
    ) where {T}
    params = generate_biot_savart_parameters()
    (; Γ,) = params
    L = params.Ls[3]  # = Lz
    filaments = generate_straight_lines(T, L, N, method)

    T_kw = BiotSavart.kelvin_wave_period(params, L)  # period of largest KWs
    dt_kw = BiotSavart.kelvin_wave_period(params, L/N)  # period of smallest resolved KWs
    dt = dt_kw

    # Forcing
    # We apply a stochastic velocity forcing.
    # The forcing changes randomly at each timestep (basic Euler–Maruyama scheme).
    # This is similar to the forcing used by Baggaley & Laurie (PRB 2014).
    rng = Random.Xoshiro()
    forcing_amplitude = Γ / 100
    forcing_ks = T[-3, -2, -1, 1, 2, 3]  # forcing wavenumbers (in the z direction)
    forcing = BrownianForcing(dt, forcing_amplitude, forcing_ks)

    forcing_velocity = generate_forcing_function(forcing)
    forcing_streamfunction = generate_streamfunction_function(forcing)

    ## Callback function to be called after each timestep
    times = T[]
    energy = T[]
    line_length = T[]

    function callback(iter)
        (; nstep, t,) = iter
        Tend = iter.prob.tspan[2]
        if nstep == 0  # make sure vectors are empty at the beginning of the simulation
            empty!(times)
            empty!(energy)
            empty!(line_length)
        end
        # Set random amplitudes and phases of forcing for next iteration.
        update_forcing!(rng, forcing)
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
    Random.seed!(rng, 42)  # restart the random number generator
    init_forcing!(rng, forcing)
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
    Random.seed!(rng, 42)  # restart the random number generator
    init_forcing!(rng, forcing)
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
