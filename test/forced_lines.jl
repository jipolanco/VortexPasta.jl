# Test initially straight infinite lines subject to velocity forcing.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm
using Random: Random
using StableRNGs: StableRNG
using FFTW: fft!, bfft!, fftfreq
using VortexPasta.PredefinedCurves: define_curve, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using UnicodePlots: lineplot, lineplot!
using JET: JET
using KernelAbstractions: KernelAbstractions as KA  # for JET only

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

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
struct DeltaCorrelatedForcing{T} <: AbstractStochasticForcing{T}
    dt :: T          # simulation timestep
    σ² :: T          # fixed forcing amplitude -- has units of a diffusion coefficient [L²T⁻¹]
    amplitude :: T   # = √{σ² / dt} -- has units of a velocity [LT⁻¹]
    ks :: Vector{T}  # forced wavenumbers
    ws :: Vector{Complex{T}}  # random forcing coefficients (~ Normal(0, std = amplitude))
end

correlation_time(f::DeltaCorrelatedForcing) = f.dt

function DeltaCorrelatedForcing(dt, σ², ks)
    ws = similar(ks, complex(eltype(ks)))
    fill!(ws, 0)
    amplitude = sqrt(σ² / dt)
    DeltaCorrelatedForcing(dt, σ², amplitude, ks, ws)
end

function update_forcing!(rng::Random.AbstractRNG, f::DeltaCorrelatedForcing)
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
        ws[i] = (1 - dt / τ) * w + amplitude * Random.randn(rng, typeof(w))
    end
    f
end

# ================================================================================ #

function dissipate_fourier!(ws, ks; dt, α, ν,)
    for i ∈ eachindex(ws, ks)
        k = ks[i]
        iszero(k) && continue  # avoid division by zero
        k² = k * k
        k⁴ = k² * k²
        ws[i] = ws[i] * exp(-dt * (α / k² + ν * k⁴))
    end
    ws
end

function dissipate!(f::ClosedFilament, L; kwargs...)
    Xs = nodes(f)
    ws = @. getindex(Xs, 1) + im * getindex(Xs, 2)
    fft!(ws)  # in-place FFT
    N = length(Xs)
    @. ws = ws / N  # normalise FFT
    ks = fftfreq(N, 2π * N / L)  # wavenumbers associated to FFT
    dissipate_fourier!(ws, ks; kwargs...)
    # Apply backwards FFT and modify `f` the filament with the results.
    bfft!(ws)
    for i ∈ eachindex(f, ws)
        f[i] = (real(ws[i]), imag(ws[i]), f[i].z)
    end
    f
end

function redistribute_along_z!(f::AbstractFilament; rtol = 1e-12, nmax = 10)
    N = length(f)
    ts = knots(f)
    z₀ = f[1].z
    Lz = f[end + 1].z - z₀
    atol = rtol * abs(Lz)  # absolute tolerance in terms of the filament period
    @assert atol > 0
    zs_target = range(z₀, z₀ + Lz; length = N + 1)
    for i ∈ eachindex(f)[2:end]
        z_target = zs_target[i]
        tn = ts[i]  # start from the node to be modified
        s⃗ = f[i]
        err = abs(s⃗.z - z_target)
        # Perform a few Newton iterations.
        # This usually converges quite fast, usually in about 2 iterations.
        n = 0
        while err > atol && n < nmax
            n += 1
            z′ = f(tn, Derivative(1)).z
            tn = tn + (z_target - s⃗.z) / z′
            s⃗ = f(tn)
            err = abs(s⃗.z - z_target)
        end
        f[i] = s⃗
    end
    f
end

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
    # Parameters are relevant for HeII if we interpret dimensions in cm and s.
    L = 2π
    Γ = 9.97e-4
    a = 1e-8
    Δ = 1/4
    Ngrid = 16
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    Ls = (L, L, L)
    β = 3.5
    α = kmax / 2β
    rcut = β / α
    ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
        # backend_short = CellListsBackend(2),
        backend_short = NaiveShortRangeBackend(),  # needed when rcut > L/2 (which is the case here, since Ngrid is small)
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4), σ = 1.5),
        quadrature = GaussLegendre(3),
    )
end

function test_forced_lines(
        ::Type{T} = Float64;
        N = 32,
        kf_norm = [1, 2, 3],  # forced wavenumbers (absolute values only)
        method = QuinticSplineMethod(),
        scheme = RK4(),
    ) where {T}
    all(>(0), kf_norm) || throw(ArgumentError("kf_norm can only contain positive values"))
    params = generate_biot_savart_parameters()
    (; Γ,) = params
    L = params.Ls[3]  # = Lz
    filaments = generate_straight_lines(T, L, N, method)

    dt_kw = BiotSavart.kelvin_wave_period(params, L/N)  # period of smallest resolved KWs
    dt = dt_kw  # simulation timestep

    # Forcing
    # We apply a stochastic velocity forcing.
    # The forcing changes randomly at each timestep (basic Euler–Maruyama scheme).
    # This is similar to the forcing used by Baggaley & Laurie (PRB 2014).
    rng = StableRNG(42)
    forcing_amplitude = Γ / 200  # this has the units of a diffusion coefficient [L²T⁻¹]
    forcing_ks = let kf = convert(Vector{T}, kf_norm)  # forcing wavenumbers (in the z direction)
        append!(kf, -kf)  # also force negative wavenumbers
    end
    forcing_τ = 4 * dt_kw
    forcing = DeltaCorrelatedForcing(dt, forcing_amplitude, forcing_ks)
    # forcing = OrnsteinUhlenbeckForcing(dt, forcing_amplitude, forcing_τ, forcing_ks)

    forcing_velocity = generate_forcing_function(forcing)
    forcing_streamfunction = generate_streamfunction_function(forcing)

    # Dissipation parameters
    kf_min = minimum(kf_norm)  # corresponds to the largest forcing scale
    kmax = π * N / L
    dissipate_α = Γ * kf_min^4 * 2  # [L⁻²T⁻¹]
    dissipate_ν = Γ / kmax^2 * 2    # [L⁴T⁻¹]

    # Callback function to be called after each timestep
    times = T[]
    energy = T[]
    line_length = T[]

    function affect!(iter)
        # Set random amplitudes and phases of forcing for next iteration.
        update_forcing!(rng, forcing)
        # Add dissipation and analyse Fourier modes
        for f ∈ iter.fs
            redistribute_along_z!(f)  # make sure nodes are uniformly distributed in z
            dissipate!(f, L; dt = iter.dt, α = dissipate_α, ν = dissipate_ν)
            update_coefficients!(f)   # update interpolation coefficients
        end
        nothing
    end

    function callback(iter)
        (; nstep, t,) = iter
        Tend = iter.prob.tspan[2]
        if nstep == 0  # make sure vectors are empty at the beginning of the simulation
            empty!(times)
            empty!(energy)
            empty!(line_length)
        end
        E = Diagnostics.kinetic_energy_from_streamfunction(iter)
        Lvort = Diagnostics.filament_length(iter)
        push!(times, t)
        push!(energy, E)
        push!(line_length, Lvort)
        # @show nstep, t/Tend, E, Lvort/4
        nothing
    end

    ## Initialise solver
    tspan = (0.0, 10 * dt)  # run a few timesteps
    prob = VortexFilamentProblem(filaments, tspan, params)

    # 1. We first test *without* the external streamfunction, which means that energy estimates
    # will only include the energy associated to the Biot–Savart velocity (which is probably
    # the most relevant energy).
    Random.seed!(rng, 42)  # restart the random number generator
    init_forcing!(rng, forcing)
    iter = init(
        prob, scheme;
        dt,
        callback, affect!,
        external_velocity = forcing_velocity,
        # external_streamfunction = forcing_streamfunction,
    )
    if VERBOSE
        @time solve!(iter)
    else
        solve!(iter)
    end
    if VERBOSE
        plt = lineplot(times, line_length ./ first(line_length); title = "Forced lines", name = "Line length", xlabel = "Time", ylabel = "L/L₀")
        lineplot!(plt, times, energy ./ first(energy); name = "Energy")
        display(plt)
    end
    L_growth = last(line_length) / first(line_length) - 1
    E_growth = last(energy) / first(energy) - 1
    VERBOSE && @show L_growth E_growth
    @test 0.005 < L_growth < 0.006
    @test 0.005 < E_growth < 0.006

    energy_prev = copy(energy)
    line_length_prev = copy(line_length)

    # 2. Same thing without dissipation (without `affect!`). In this case line length and
    # energy should be larger than before.
    Random.seed!(rng, 42)  # restart the random number generator
    init_forcing!(rng, forcing)
    iter = init(
        prob, scheme;
        dt,
        callback,
        # affect!,  # commented => no dissipation
        external_velocity = forcing_velocity,
        # external_streamfunction = forcing_streamfunction,
    )
    if VERBOSE
        @time solve!(iter)
    else
        solve!(iter)
    end
    L_growth = last(line_length) / first(line_length) - 1
    E_growth = last(energy) / first(energy) - 1
    VERBOSE && @show L_growth E_growth
    @test 0.035 < L_growth < 0.045
    @test 0.030 < E_growth < 0.040

    # 3. Now test including the external_streamfunction parameter (energy values will be much
    # higher). Note that this doesn't change the dynamics, and the vortex trajectories
    # should be exactly the same as before.
    Random.seed!(rng, 42)  # restart the random number generator
    init_forcing!(rng, forcing)
    iter = init(
        prob, scheme;
        dt,
        callback, affect!,
        external_velocity = forcing_velocity,
        external_streamfunction = forcing_streamfunction,
    )
    if VERBOSE
        @time solve!(iter)
    else
        solve!(iter)
    end
    @test line_length ≈ line_length_prev  # we get the exact same results
    # Energy estimates are much larger when including energy associated to the
    # forcing velocity.
    # @show energy ./ energy_prev
    @test all(8 .< energy ./ energy_prev .< 40)

    ## Inference tests
    enable_jet = get(ENV, "JULIA_ENABLE_JET_KA_TESTS", "false") ∈ ("true", "1")  # enable JET tests involving KA kernels
    if enable_jet
        # Check that everything is inferred.
        # (This is only for tests and can be ignored!)
        JET.@test_opt ignored_modules=(Base, Base.PCRE, KA, Base.IteratorsMD) init(
            prob, scheme;
            dt,
            callback, affect!,
            external_velocity = forcing_velocity,
            external_streamfunction = forcing_streamfunction,
        )
        JET.@test_opt ignored_modules=(Base, Base.PCRE, KA, Base.IteratorsMD) step!(iter)
    end

    nothing
end

@testset "Velocity forcing" begin
    test_forced_lines()
end
