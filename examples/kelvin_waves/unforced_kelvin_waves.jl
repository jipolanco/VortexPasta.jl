# Simulate nearly straight infinite lines subject to a small initial perturbation.

using Random: Random
using FFTW: fft!, bfft!, fftfreq
using StaticArrays: SVector
using VortexPasta.PredefinedCurves: define_curve, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using UnicodePlots: lineplot, lineplot!


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

function generate_straight_lines(::Type{T}, L, N, method; p::PeriodicLine = PeriodicLine()) where {T}
    funcs = let L = L
        (
            define_curve(p; scale = (+1, +1, +L), translate = (1L/4, 1L/4, L/2), orientation = +1),
            define_curve(p; scale = (-1, +1, +L), translate = (3L/4, 1L/4, L/2), orientation = -1),
            define_curve(p; scale = (+1, -1, +L), translate = (1L/4, 3L/4, L/2), orientation = -1),
            define_curve(p; scale = (-1, -1, +L), translate = (3L/4, 3L/4, L/2), orientation = +1),
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
    Ngrid = floor(Int, 32 * 4/5)
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    Ls = (L, L, L)
    α = kmax / 5
    rcut = 5 / α
    ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
        # backend_short = CellListsBackend(2),
        backend_short = NaiveShortRangeBackend(),  # needed when rcut > L/2 (which is the case here, since Ngrid is small)
        # backend_long = NonuniformFFTsBackend(; σ = 2.0, m = HalfSupport(4)),
        backend_long = NonuniformFFTsBackend(),
        # backend_long = FINUFFTBackend(),
        quadrature = GaussLegendre(2),
    )
end

function run_unforced_lines(
        ::Type{T} = Float64;
        N = 128,
        tmax = 0.5,  # maximum time in units of T_kw
        A = 0.01,
        ks_init = (1,),
        method = FourierMethod(),
        scheme = RK4(),
    ) where {T}
    params = generate_biot_savart_parameters()
    L = params.Ls[3]  # = Lz
    println(params)

    # Generate perturbed lines
    rng = Random.Xoshiro(42)
    filaments = let ks = SVector(T.(ks_init))
        phases = 2 * rand(rng, typeof(ks))  # in [0, 2]
        amplitudes = randn(rng, typeof(ks))
        amplitudes = amplitudes * A^2 / sum(abs2, amplitudes)  # renormalise
        @info "Perturbation parameters:"
        @show ks phases amplitudes
        r_pert(t) = sum(eachindex(ks)) do i
            amplitudes[i] * cospi(2 * ks[i] * t + phases[i])
        end
        generate_straight_lines(T, L, N, method; p = PeriodicLine(r = r_pert))
    end

    # Simulation timestep
    T_kw = BiotSavart.kelvin_wave_period(params, L)  # period of largest KWs
    dt_kw = BiotSavart.kelvin_wave_period(params, L/N)  # period of smallest resolved KWs
    dt = dt_kw / 2

    # Callback function to be called after each timestep
    times = T[]
    energy = T[]
    line_length = T[]
    ks_pos = T[]

    ws_time = Complex{T}[]

    function affect!(iter)
        for f ∈ iter.fs
            redistribute_along_z!(f)  # make sure nodes are uniformly distributed in z
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
        let f = iter.fs[2]
            local N = length(f)
            sizehint!(ws_time, length(ws_time) + N)
            for s⃗ ∈ nodes(f)
                push!(ws_time, s⃗.x + im * s⃗.y)
            end
        end
        if nstep % 10 == 0
            @show nstep, t/Tend, E, L_growth
            # write_vtkhdf("vtk/kws_$nstep.vtkhdf", iter.fs) do vtk
            #     vtk["velocity"] = iter.vs
            # end
        end
        nothing
    end

    # Initialise solver
    Tend = tmax * T_kw
    @show T_kw, Tend, dt, Tend/dt
    tspan = (0.0, Tend)
    prob = VortexFilamentProblem(filaments, tspan, params)

    iter = init(
        prob, scheme;
        dt,
        callback, affect!,
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
    @show L_growth E_growth  # -8e-14, -7e-14

    # Return all results
    results = (;
        N, times, params, energy, line_length,
        ks_pos, ws_time,
    )

    results
end

##

tmax = let params = generate_biot_savart_parameters(), kinit = 10
    T_kw = BiotSavart.kelvin_wave_period(params, 2π)
    t_k = BiotSavart.kelvin_wave_period(params, 2π / kinit)
    40 * t_k / T_kw
end

results = run_unforced_lines(;
    N = 64, tmax, A = 2π / 1000,
    ks_init = (4, 5, 6),
    # ks_init = (10,),
    method = FourierMethod(),
);

##

# Plot spatiotemporal spectrum of filament fluctuations

using GLMakie
using FFTW: fft!, fft, bfft, fftfreq, fftshift

(; times, N,) = results
inds = eachindex(times)[end÷4:end]  # drop the initial 1/4 of the simulation
ws = reshape(results.ws_time, N, :)[:, inds]
Nt = length(times[inds])
Nt_out = size(ws, 2)
step_out = ceil(Int, Nt / Nt_out)
dt_out = times[begin + step_out] - times[begin]
ks = fftfreq(N, N)
ωs = fftfreq(Nt_out, 2π / dt_out)
ws_h = copy(ws)
ws_mean = sum(i -> ws[i, 1], 1:N) / N
ws_h .-= ws_mean
fft!(ws_h)
ws_h ./= length(ws_h)

amplitudes = fftshift(log10.(abs2.(ws_h)))
ks_shift = fftshift(ks)
ωs_shift = fftshift(ωs)

ks_pos = ks[2:(end ÷ 2)]
ωs_kw = 2π ./ BiotSavart.kelvin_wave_period.(Ref(results.params), 2π ./ ks_pos)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "k", ylabel = "ω")
hm = heatmap!(
    ax, ks_shift, ωs_shift, amplitudes;
    colorrange = (-34, -12),
)
lines!(ax, ks_pos, ωs_kw; color = (:white, 0.5), linestyle = :dash)
Colorbar(fig[1, 2], hm)
DataInspector(fig)
fig

##

# Plot average spatial spectrum of filament fluctuations

nk = zeros(length(ks_pos))
# ws_t = bfft(ws_h, 2)
ws_t = fft(ws, 1) ./ N
let n = 0
    for j ∈ axes(ws_t, 2)
        for i ∈ eachindex(nk)
            nk[i] += abs2(ws_t[begin + i, j]) + abs2(ws_t[end - i + 1, j])
        end
        n += 1
    end
    nk ./= n
end

fig = Figure()
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10)
scatterlines!(ax, ks_pos, nk)
let ks = 10.0.^(range(log10(5), log10(20); length = 3))
    lines!(ax, ks, 1e-20 * ks.^(-11/3); color = :grey, linestyle = :dash)
end
fig
