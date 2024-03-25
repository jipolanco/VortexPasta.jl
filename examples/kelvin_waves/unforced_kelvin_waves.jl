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

function generate_biot_savart_parameters(::Type{T}) where {T}
    # Parameters are relevant for HeII if we interpret dimensions in cm and s.
    cm = 1.0    # one cm = 1 length unit
    sec = 1.0   # one second = 1 time unit
    L = 2π      # period in length units
    Γ = 9.97e-4 * cm^2/sec
    a = 1e-8 * cm
    Δ = 1/4
    Ngrid = floor(Int, 32)
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    Ls = (L, L, L)
    β = 4.5
    α = kmax / 2β
    rcut = β / α
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        # backend_short = CellListsBackend(2),
        backend_short = NaiveShortRangeBackend(),  # needed when rcut ≈ L/2 (which is the case here, since Ngrid is small)
        backend_long = NonuniformFFTsBackend(; σ = 1.5, m = HalfSupport(8)),
        # backend_long = FINUFFTBackend(tol = 1e-8),
        # backend_long = ExactSumBackend(),
        quadrature = GaussLegendre(4),
    )
end

function run_unforced_lines(
        params::ParamsBiotSavart;
        N = 128,
        tend = 10.0,  # final time
        A = 0.01,
        ks_init::AbstractVector = [1],  # perturbed wavenumbers (absolute values only)
        method = FourierMethod(),
        scheme = RK4(),
        dt_factor = 1.0,
    )
    all(>(0), ks_init) || throw(ArgumentError("ks_init can only contain positive values"))
    T = eltype(params)  # numerical precision (e.g. Float64)
    L = params.Ls[3]  # = Lz
    println(params)

    # Generate perturbed lines
    rng = Random.Xoshiro(42)
    ks_init_all = let kf = convert(Vector{T}, ks_init)  # initialised wavenumbers (in the z direction)
        append!(kf, -kf)  # also initialise negative wavenumbers
        kf
    end
    filaments = let ks = ks_init_all
        # α = -2
        # α = -11/3
        α = -5
        ws = randn(rng, Complex{T}, length(ks))  # random complex perturbation (includes random amplitudes + phases)
        for i ∈ eachindex(ks, ws)
            ws[i] *= abs(ks[i])^(α/2)
        end
        ws .*= A / sqrt(sum(abs2, ws))           # renormalise to get the wanted perturbation amplitude
        @assert sum(abs2, ws) ≈ A^2
        r_pert(t) = sum(eachindex(ks)) do i
            ws[i] * cispi(2 * ks[i] * t)
        end
        generate_straight_lines(T, L, N, method; p = PeriodicLine(r = r_pert))
    end

    # Simulation timestep
    dt_kw = BiotSavart.kelvin_wave_period(params, L/N)  # period of smallest resolved KWs
    dt = dt_kw * dt_factor

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
        if nstep % 100 == 0
            @show nstep, t/Tend, E, L_growth
            # write_vtkhdf("vtk/kws_$nstep.vtkhdf", iter.fs) do vtk
            #     vtk["velocity"] = iter.vs
            # end
        end
        nothing
    end

    # Initialise solver
    @show tend, dt, tend/dt
    tspan = (0.0, tend)
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
        iter,
        N, times, params, energy, line_length,
        ks_pos, ws_time,
    )

    results
end

##

params = generate_biot_savart_parameters(Float64)

# This roughly corresponds to the initial perturbation used by Krstulovic (PRE 2012) in GP
# simulations (runs II, III).
# In the paper, the first two modes are forced. The perturbation amplitude is A = 4ξ (with
# ξ the healing length), and ξ = Δx = L/N⟂ ≈ 0.025, where L = 2π, and N⟂ = 256 is the resolution in
# the perpendicular directions. So, in terms of the wavelength λ_max of the largest initialised
# KW, the amplitude is A = 4 * λ_max / 256.
ks_init = [1, 2]
Lz = params.Ls[3]
λ_max = Lz / maximum(ks_init)
A = 4 * λ_max / 256

T_kw = BiotSavart.kelvin_wave_period(params, λ_max)  # period of slowest Kelvin waves
tend = 1.0 * T_kw  # simulation time

results = run_unforced_lines(
    params;
    N = 128, tend, A,
    ks_init,
    method = FourierMethod(),
    # scheme = Strang(RK4(); nsubsteps = 1),
    scheme = RK4(),
    dt_factor = 1.0,
);

## Plot spatiotemporal spectrum of filament fluctuations (dispersion relation)

using GLMakie
using FFTW: fft!, fft, bfft, fftfreq, fftshift
using DSP: DSP

(; times, N, ws_time,) = results
inds = eachindex(times)[2end÷4:end]  # drop the initial half of the simulation
ws = reshape(ws_time, N, :)[:, inds]
ts = times[inds] ./ T_kw  # normalised time
Nt = length(times[inds])
Nt_out = size(ws, 2)
step_out = ceil(Int, Nt / Nt_out)
dt_out = times[begin + step_out] - times[begin]
ks = fftfreq(N, N)
ωs = fftfreq(Nt_out, 2π / dt_out)
ws_h = copy(ws)
ws_mean = sum(i -> ws[i, 1], 1:N) / N
ws_h .-= ws_mean

# Apply window function in the temporal axis to make the signal closed to being time-periodic.
# See Kelvin wave tutorial for more details (in 1D).
window = DSP.Windows.hanning(size(ws_h, 2))
for i ∈ axes(ws_h, 1)
    @. @views ws_h[i, :] *= window
end

fft!(ws_h)
ws_h ./= length(ws_h)

amplitudes = fftshift(log10.(abs2.(ws_h)))
ks_shift = fftshift(ks)
ωs_shift = fftshift(ωs)

ks_pos = ks[2:(end ÷ 2)]
ωs_kw = 2π ./ BiotSavart.kelvin_wave_period.(Ref(results.params), 2π ./ ks_pos)

Makie.set_theme!(
    Axis = (
        xlabelsize = 20, ylabelsize = 20,
        xticklabelsize = 16, yticklabelsize = 16,
    ),
)

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"k", ylabel = L"ω")
hm = heatmap!(
    ax, ks_shift, ωs_shift, amplitudes;
    colorrange = (-40, -10),
)
lines!(ax, ks_pos, ωs_kw; color = (:white, 0.2), linestyle = :dash)
Colorbar(fig[1, 2], hm; label = L"\log \, |\hat{w}(k,\,ω)|^2", labelsize = 20)
DataInspector(fig)
fig

## Plot spectra over time

m = length(times) ÷ 20  # plot ~20 times
m = max(2, m)
ts = times ./ T_kw  # normalised time
inds = eachindex(ts)[m:m:end]
ts_lims = (ts[begin], ts[end])
ws = reshape(ws_time, N, :)

nk_avg = zeros(length(ks_pos))  # average spectrum
nk = similar(nk_avg)
cmap = Makie.to_colormap(:roma)

fig = Figure(size = (1200, 500))
ax_normal = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = L"k", ylabel = L"n(k)")
ax_compen = Axis(fig[1, 2]; xscale = log10, yscale = log10, xlabel = L"k", ylabel = L"k^{11/3} \, n(k)")
for i ∈ inds
    local wi = @. ws[:, i]
    local ŵi = fft(wi) ./ N
    for i ∈ eachindex(nk)
        nk[i] = abs2(ŵi[begin + i]) + abs2(ŵi[end - i + 1])
    end
    nk_avg .+= nk
    t = ts[i]
    color = Makie.interpolated_getindex(cmap, t, ts_lims)
    label = string(round(t; digits = 2))
    scatterlines!(ax_normal, ks_pos, nk; color, label)
    let ys = @. ks_pos^(11/3) * nk
        scatterlines!(ax_compen, ks_pos, ys; color, label)
    end
end
nk_avg ./= length(inds)
let kws = (color = :orangered, linewidth = 3, label = "Mean")
    lines!(ax_normal, ks_pos, nk_avg; kws...)
    let ys = @. ks_pos^(11/3) * nk_avg
        scatterlines!(ax_compen, ks_pos, ys; kws...)
    end
end
let ks = 10.0.^(range(log10(last(ks_init) * 2), log10(30); length = 3))
    lines!(ax_normal, ks, 1e-13 * ks.^(-11/3); color = :grey, linestyle = :dash)
end
cb = Colorbar(fig[1, end + 1]; colormap = cmap, colorrange = ts_lims, label = L"Time $(t / T_1)$", labelsize = 20)
DataInspector(fig)
fig

