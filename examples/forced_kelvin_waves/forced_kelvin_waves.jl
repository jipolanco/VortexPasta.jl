# Simulate initially straight infinite lines subject to velocity forcing.

using Random: Random
using FFTW: fft!, bfft!, fftfreq
using VortexPasta.PredefinedCurves: define_curve, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using UnicodePlots: lineplot, lineplot!
using GLMakie

# ================================================================================ #
## Forcing definitions

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

## Delta-correlated forcing (like Brownian motion)
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

## Stochastic forcing with correlation time τ.
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
## Fourier space stuff

function wave_action_spectrum(ks::AbstractVector, rhat::AbstractVector)
    @assert ks[2] == -ks[end]  # contains positive and negative wavenumbers
    @assert length(ks) == length(rhat)
    N = length(ks)
    if iseven(N)
        Nh = N ÷ 2
        @assert ks[Nh + 1] ≈ -(ks[Nh] + ks[2])  # wavenumbers change sign after index Nh
    else
        Nh = N ÷ 2 + 1
        @assert ks[Nh + 1] == -ks[Nh]  # wavenumbers change sign after index Nh
    end
    ks_pos = ks[2:Nh]  # only positive wavenumbers
    nk = similar(ks_pos)
    for j ∈ eachindex(ks_pos)
        local k = ks_pos[j]
        i⁺ = 1 + j      # index of coefficient corresponding to wavenumber +k
        i⁻ = N + 1 - j  # index of coefficient corresponding to wavenumber -k
        @assert ks[i⁺] == -ks[i⁻] == k  # verification
        nk[j] = abs2(rhat[i⁺]) + abs2(rhat[i⁻])
    end
    ks_pos, nk
end

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

function dissipate_and_analyse!(f::ClosedFilament, L; kwargs...)
    Xs = nodes(f)
    ws = @. getindex(Xs, 1) + im * getindex(Xs, 2)
    fft!(ws)  # in-place FFT
    N = length(Xs)
    @. ws = ws / N  # normalise FFT
    ks = fftfreq(N, 2π * N / L)  # wavenumbers associated to FFT
    dissipate_fourier!(ws, ks; kwargs...)
    ks_pos, nk = wave_action_spectrum(ks, ws)
    # Apply backwards FFT and modify `f` the filament with the results.
    bfft!(ws)
    for i ∈ eachindex(f, ws)
        f[i] = (real(ws[i]), imag(ws[i]), f[i].z)
    end
    (; ks_pos, nk,)
end

# ================================================================================ #

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
    α = kmax / 5
    rcut = 5 / α
    ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
        # backend_short = CellListsBackend(2),
        backend_short = NaiveShortRangeBackend(),  # needed when rcut > L/2 (which is the case here, since Ngrid is small)
        backend_long = FINUFFTBackend(),
        quadrature = GaussLegendre(2),
    )
end

function run_forced_lines(
        ::Type{T} = Float64;
        N = 128,
        tmax = 0.5,        # maximum time in units of T_kw
        kf_norm = [2, 3],  # forced wavenumbers (absolute values only)
        method = QuinticSplineMethod(),
        scheme = RK4(),
    ) where {T}
    all(>(0), kf_norm) || throw(ArgumentError("kf_norm can only contain positive values"))
    params = generate_biot_savart_parameters()
    (; Γ,) = params
    L = params.Ls[3]  # = Lz
    filaments = generate_straight_lines(T, L, N, method)
    println(params)

    # Simulation timestep (also appearing in the amplitude of the stochastic forcing)
    T_kw = BiotSavart.kelvin_wave_period(params, L)  # period of largest KWs
    dt_kw = BiotSavart.kelvin_wave_period(params, L/N)  # period of smallest resolved KWs
    dt = dt_kw

    # Forcing
    # We apply a stochastic velocity forcing.
    # The forcing changes randomly at each timestep (basic Euler–Maruyama scheme).
    # This is similar to the forcing used by Baggaley & Laurie (PRB 2014).
    rng = Random.Xoshiro()
    forcing_amplitude = Γ * 1e-2  # this has the units of a diffusion coefficient [L²T⁻¹]
    forcing_ks = let kf = convert(Vector{T}, kf_norm)  # forcing wavenumbers (in the z direction)
        append!(kf, -kf)  # also force negative wavenumbers
    end
    forcing_τ = 4 * dt_kw
    @show forcing_amplitude forcing_ks forcing_τ
    @assert forcing_τ > dt
    # forcing = DeltaCorrelatedForcing(dt, forcing_amplitude, forcing_ks)
    forcing = OrnsteinUhlenbeckForcing(dt, forcing_amplitude, forcing_τ, forcing_ks)

    # Dissipation parameters
    kf_min = minimum(kf_norm)  # corresponds to the largest forced scale
    kmax = π * N / L
    dissipate_α = Γ * kf_min^4 * 2  # [L⁻²T⁻¹]
    dissipate_ν = Γ / kmax^2 * 2   # [L⁴T⁻¹]
    @show dissipate_α dissipate_ν

    # Callback function to be called after each timestep
    times = T[]
    energy = T[]
    line_length = T[]
    nk_time = Vector{T}[]  # vector of vectors
    forcing_coef = Complex{T}[]  # just to visualise the time evolution of the forcing coefficients
    ks_pos = T[]

    function affect!(iter)
        # Set random amplitudes and phases of forcing for next iteration.
        update_forcing!(rng, forcing)
        push!(forcing_coef, forcing.ws[1])

        # Add dissipation and analyse Fourier modes
        for (i, f) ∈ enumerate(iter.fs)
            redistribute_along_z!(f)  # make sure nodes are uniformly distributed in z
            data = dissipate_and_analyse!(f, L; dt = iter.dt, α = dissipate_α, ν = dissipate_ν)
            update_coefficients!(f; knots = knots(f))  # update interpolation coefficients
            if i == 1  # only keep data for a single filament
                push!(nk_time, data.nk)
                copy!(ks_pos, data.ks_pos)
            end
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
        local L_growth = Lvort/4L - 1  # relative growth of vortex lines
        if nstep % 10 == 0
            @show nstep, t/Tend, E, L_growth
            write_vtkhdf("vtk/kws_$nstep.hdf", iter.fs) do vtk
                vtk["velocity"] = iter.vs
            end
        end
        nothing
    end

    ## Initialise solver
    Tend = tmax * T_kw
    @show T_kw, correlation_time(forcing), Tend, dt, Tend/dt
    tspan = (0.0, Tend)
    prob = VortexFilamentProblem(filaments, tspan, params)

    forcing_velocity = generate_forcing_function(forcing)

    Random.seed!(rng, 42)
    init_forcing!(rng, forcing)
    iter = init(
        prob, scheme;
        dt,
        callback, affect!,
        external_velocity = forcing_velocity,
    )
    @time solve!(iter)
    println(iter.to)
    let
        plt = lineplot(times, line_length ./ first(line_length) .- 1; title = "Forced lines", name = "Line length", xlabel = "Time", ylabel = "L/L₀ - 1")
        lineplot!(plt, times, energy ./ first(energy) .- 1; name = "Energy")
        display(plt)
    end
    L_growth = last(line_length) / first(line_length) - 1
    E_growth = last(energy) / first(energy) - 1
    @show L_growth E_growth

    let
        plt = lineplot(times, real.(forcing_coef))
        lineplot!(plt, times, imag.(forcing_coef))
        display(plt)
    end

    # Return all results
    results = (;
        times, energy, line_length,
        ks_pos, nk_time, forcing_coef,
    )

    results
end

##

mkpath("vtk")  # for VTK output
results = run_forced_lines(; N = 128, tmax = 0.1,)

## Plot results

let fig = Figure()
    ax = Axis(fig[1, 1]; xscale = log10, yscale = log10,)
    step = 20
    inds = eachindex(results.times)[begin:step:end]
    tmax = results.times[end]
    for i ∈ inds
        ks = results.ks_pos
        t = results.times[i]
        nk = results.nk_time[i]
        all(iszero, nk) && continue  # skip if all values are zero (typically at t = 0)
        cmap = Makie.to_colormap(:viridis)
        color = Makie.interpolated_getindex(cmap, t/tmax)  # colour as a function of time
        scatterlines!(ax, ks, nk; label = "$t", color,)
    end
    let ks = 10.0.^range(log10(2), log10(16); length = 3)
        lines!(ax, ks, ks.^(-11/3); color = :grey, linestyle = :dash)
    end
    DataInspector(fig)
    fig
end
