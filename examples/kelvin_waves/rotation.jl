# Verify the Kelvin wave dispersion relation for a superfluid under rotation (including
# collective effects).

using VortexPasta
using VortexPasta.BiotSavart
using VortexPasta.Diagnostics
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.FilamentIO
using VortexPasta.Quadratures: Quadratures, GaussLegendre
using VortexPasta.PredefinedCurves: PeriodicLine, define_curve
using VortexPasta.Timestepping
using StableRNGs: StableRNG  # for deterministic random number generation
using Bessels: besselk0  # used in dispersion relation of KWs under rotation
using FFTW: fft, bfft, fft!, bfft!, fftfreq, fftshift
using DSP: DSP  # for Fourier analysis / window functions
using Statistics: std
using GLMakie, CairoMakie, MathTeXEngine  # for plots
using InverseFunctions: InverseFunctions

# Set TeX Gyre Pagella as text and math font.
MathTeXEngine.set_texfont_family!(FontFamily("TeXGyrePagella"))
Makie.set_theme!(theme_latexfonts())
GLMakie.activate!()

## Helper functions

# Give a colour to a filament based on its local orientation wrt Z.
function filament_colour(f::AbstractFilament, refinement)
    cs = Float32[]
    ζs = range(0, 1; length = refinement + 1)[1:refinement]  # for interpolation
    for seg ∈ segments(f), ζ ∈ ζs
        colour = seg(ζ, UnitTangent())[3]  # in [-1, 1]
        push!(cs, colour)
    end
    let seg = last(segments(f))  # "close" the curve
        colour = seg(1.0, UnitTangent())[3]
        push!(cs, colour)
    end
    cs
end

# Plot a list of filaments
function plot_filaments(fs::AbstractVector{<:AbstractFilament}, Ls::NTuple{3})
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    rect = Rect(0, 0, 0, Ls...)
    wireframe!(ax, rect; color = (:black, 0.5), linewidth = 0.2)
    hidespines!(ax)
    for f ∈ fs
        refinement = 4
        color = filament_colour(f, refinement)
        plot!(
            ax, f;
            refinement, color, colormap = :RdBu_9, colorrange = (-1, 1), markersize = 4,
        )
    end
    fig
end

# Plot a single filament
plot_filaments(f::AbstractFilament, args...) = plot_filaments([f], args...)

# Create a random perturbation of a vortex.
# This will fill the vector `ws` with complex values corresponding to random
# perturbations in the x and y directions.
function random_perturbation!(rng, ws::AbstractVector{<:Complex}; Lz::Real, Az_rms::Real)
    N = length(ws)

    # First define perturbation in Fourier space
    ks = fftfreq(N, 2π * N / Lz)  # wavenumbers
    kmax = π * N / Lz
    # kmax_perturb = 0.5 * kmax   # perturb up to this wavenumber (to avoid issues at the discretisation scale)
    kmax_perturb = kmax   # perturb up to this wavenumber
    for i in eachindex(ws, ks)
        kabs = abs(ks[i])
        if 0 < kabs <= kmax_perturb
            ws[i] = randn(rng, ComplexF64)  # perturb all wavenumbers equally
        else
            ws[i] = 0
        end
    end

    # Rescale perturbation to obtain the wanted amplitude (rms value)
    factor = Az_rms / sqrt(sum(abs2, ws))
    ws .*= factor

    # Now transform to physical space and return the result
    bfft!(ws)  # backwards FFT

    ws
end

# Interpolate (x, y) locations on the vortex at equidistant z locations.
# Locations are written in complex form, w = x + iy.
function interpolate_at_equidistant_zs!(
        positions::AbstractVector{<:Complex}, f::AbstractFilament;
        Lz, atol = eps(Lz) * 100, maxiter = 10,
    )
    N = length(positions)
    t = first(knots(f))  # initial location on the vortex; this is usually 0
    for i in 1:N
        # Interpolate points so that they're regularly spaced (constant Δz)
        z_wanted = (i - 1) * Lz / N  # since the parameter t is equal to z
        s⃗ = f(t)
        niter = 0
        while abs(z_wanted - s⃗.z) > atol && niter < maxiter
            niter += 1
            z′ = f(t, Derivative(1)).z
            t = t + (z_wanted - s⃗.z) / z′
            s⃗ = f(t)
        end
        @assert niter < maxiter  # convergence is usually very fast
        @assert isapprox(s⃗.z, z_wanted; atol)
        w = s⃗.x + im * s⃗.y
        positions[i] = w
    end
    positions
end

# Incrementally integrate function `f` from xs[1] to xs[i] for each i.
# We use functions from the VortexPasta.Quadratures module.
function integrate_incremental!(f::F, ys, xs; quad) where {F}
    ys[begin] = 0  # integral from xs[1] to xs[1]
    @inbounds for i in eachindex(xs, ys)[2:end]
        # Integrate in (a, b) using given quadrature rule (typically GaussLegendre(n)).
        a = xs[i - 1]
        b = xs[i]
        dy = Quadratures.integrate(f, quad, (a, b))
        ys[i] = ys[i - 1] + dy
    end
    ys
end

# Perform "imaginary-time" evolution to find a stable equilibrium (minimal energy) configuration.
function find_stable_equilibrium!(xs::AbstractVector{T}, ys::AbstractVector{T}; Ls, rtol = T(1e-6), verbose = false, checkpoints = false) where {T}
    # We perform a pseudo-2D simulation with very few points in z.
    @info "Performing imaginary-time evolution"
    N = 4
    Lx, Ly, _ = Ls
    ℓ = sqrt(Lx * Ly / length(xs))  # mean inter-vortex distance
    Lz = min(Lx, Ly)  # this is not physically relevant here
    Ls = (Lx, Ly, Lz)
    ε = ℓ / 100  # random perturbation strength
    rng = StableRNG(42)
    fs = map(eachindex(xs, ys)) do i
        x, y = xs[i], ys[i]
        x += ε * randn(rng, T)
        y += ε * randn(rng, T)
        # Add small random perturbation just in case we're at an unstable equilibrium (maximal energy)
        Filaments.init(ClosedFilament{T}, N, FourierMethod()) do τ
            z = τ * Lz  # note: τ ∈ [0, 1]
            Vec3(x, y, z)
        end
    end
    rcut = min(Ls...) / 3
    Γ = one(T)
    params = ParamsBiotSavart(;
        Γ, a = T(1e-7), Δ = T(1/2),
        splitting = KaiserBesselSplitting(; Ls, β = 14.0, rcut),
        # backend_short = NaiveShortRangeBackend(),
        backend_short = CellListsBackend(),
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4)),
        quadrature = GaussLegendre(3),
    )
    if verbose
        println(params)
    end
    tsim = one(T)  # ignored
    prob = VortexFilamentProblem(fs, tsim, params)
    println(prob)
    scheme = Strang(RK4(), RK4(); nsubsteps = 4)
    # scheme = RK4()
    iter = init(prob, scheme; dt = one(T), adaptivity = AdaptBasedOnSegmentLength(2.0) | AdaptBasedOnVelocity(0.1), mode = MinimalEnergy())
    atol = rtol * Γ / ℓ  # maximum allowed velocity
    atol² = atol^2
    max_squared_velocity(iter) = maximum(vs -> maximum(v⃑ -> sum(abs2, v⃑), vs), iter.vs)
    v²_max = max_squared_velocity(iter)
    E = Diagnostics.kinetic_energy(iter)
    E_initial = E
    if verbose
        @show iter.nstep, iter.t, iter.dt, E, v²_max/atol²
    end
    tsf = TimeSeriesFile()
    if checkpoints
        let fname = "imaginary_time_step$(iter.nstep).vtkhdf"
            save_checkpoint(fname, iter)
            tsf[iter.t] = fname
        end
    end
    function check_straight(f; atol)
        x, y, _ = f[begin]
        ok = true
        for s⃗ in f
            ok = ok && isapprox(s⃗.x, x; atol)
            ok = ok && isapprox(s⃗.y, y; atol)
        end
        ok
    end
    stats_file = "imaginary_time_stats.dat"
    open(stats_file, "w") do io
        println(io, "# (1) Time  (2) Energy  (3) vL_max")
        println(io, iter.t, '\t', E, '\t', sqrt(v²_max))
    end
    while v²_max > atol²
        step!(iter)
        if !all(f -> check_straight(f; atol = ℓ * rtol), iter.fs)
            error("vortices did not remain straight during imaginary time evolution. Try reducing point spacing δ (by changing Lz or N).")
        end
        v²_max = max_squared_velocity(iter)
        E = Diagnostics.kinetic_energy(iter)
        if iter.nstep % 10 == 0
            open(stats_file, "a") do io
                println(io, iter.t, '\t', E, '\t', sqrt(v²_max))
            end
        end
        if iter.nstep % 100 == 0
            if verbose
                @show iter.nstep, iter.t, iter.dt, E, v²_max/atol²
            end
            if checkpoints
                let fname = "imaginary_time_step$(iter.nstep).vtkhdf"
                    save_checkpoint(fname, iter)
                    tsf[iter.t] = fname
                    save("imaginary_time.vtkhdf.series", tsf)
                end
            end
        end
    end
    if verbose
        println(iter.to)
    end
    E = Diagnostics.kinetic_energy(iter; quad = params.quad)
    @assert E < E_initial  # we minimised the energy
    # Now copy back the equilibrium positions.
    for i in eachindex(xs, ys, iter.fs)
        f = iter.fs[i]
        x, y, _ = f[begin]
        xs[i] = x
        ys[i] = y
    end
    (xs, ys)
end

function run_simulation(;
        nvortices = 4,       # number of vortices (typically multiple of 4)
        Ls = (2π, 2π, 2π),   # domain period
        N = 64,              # number of points per vortex (= line resolution)
        Γ = 1.0, a = 1e-7, Δ = 0,
        in_phase_perturbation = true,  # perturbations of all vortices are in phase?
        Az_rms = 1e-5 * Lz,   # amplitude of random perturbation varying in Z (rms value)
        Ax = 0 * Lx,      # amplitude of perturbation varying in X
        Ay = Ax,          # amplitude of perturbation varying in Y
        kx = 0, ky = 0,       # mode of horizontal perturbations
        Axy_rms = 0 * Lx,  # amplitude of random perturbation in the horizontal directions
        method = QuinticSplineMethod(),  # filament discretisation method
        rcut_factor = 2/5,
        dt_factor = 1.0,
        nrotations = 4,  # number of rotation periods to simulate
        nsave = 100,     # number of checkpoints to save
    )
    #== Biot-Savart parameters ==#
    rcut = min(Ls...) * rcut_factor
    T = typeof(Γ)  # usually Float64
    params = ParamsBiotSavart(
        T;
        Γ, a, Δ,
        # splitting = GaussianSplitting(; Ls, β = 3.5, rcut),
        splitting = KaiserBesselSplitting(; Ls, β = 14.0, rcut),
        quadrature = GaussLegendre(3),        # quadrature for integrals over filament segments
        quadrature_near_singularity = GaussLegendre(3),
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4)),
        backend_short = CellListsBackend(2),
    )
    println('\n', params)

    #== Initialise vortices ==#
    # Create "base" vortex from which we will create the actual vortices.
    # This is a straight vortex that passes through the origin and goes from z = 0 to z = Lz.
    p = PeriodicLine()
    S = define_curve(p; scale = (1, 1, Lz), translate = (0, 0, 0.5 * Lz))
    ts = range(0, 1; length = N + 1)[1:N]  # evaluation points in [0, 1) -- we make sure that z = 0 is included
    f_base = Filaments.init(S, ClosedFilament, ts, method)
    # plot_filaments(f_base, params.Ls)

    # Create actual vortices, including random perturbation
    fs = [f_base]; empty!(fs)  # creates empty vector with the right element type
    ws = zeros(ComplexF64, N)  # this will be used to store random perturbations
    rng = StableRNG(42)        # initialise random number generator (RNG)

    # Determine stable equilibrium positions using imaginary-time evolution.
    # We start from random vortex positions.
    xs_base = Ls[1] * rand(rng, nvortices)
    ys_base = Ls[2] * rand(rng, nvortices)
    find_stable_equilibrium!(xs_base, ys_base; Ls, verbose = true, checkpoints = true)
    w_equilibrium = complex.(xs_base, ys_base)  # equilibrium position (x + iy) of each filament

    # Apply random perturbations in XY plane
    if Axy_rms > 0
        xs_base .+= Axy_rms * randn(rng, nvortices)
        ys_base .+= Axy_rms * randn(rng, nvortices)
    end

    for i in eachindex(xs_base, ys_base)
        f = similar(f_base)
        if !in_phase_perturbation || isempty(fs)
            # Don't update ws if in_phase_perturbation == true and this is not the first vortex.
            random_perturbation!(rng, ws; Lz, Az_rms)  # compute random perturbation, modifying ws
        end
        x₀ = xs_base[i]  # x position of vortex (before perturbation in Z)
        y₀ = ys_base[i]  # y position of vortex (before perturbation in Z)
        dx = Ax * sinpi(2 * kx * x₀ / Lx)
        dy = Ay * sinpi(2 * ky * y₀ / Ly)
        for n in eachindex(f, f_base)
            x = x₀ + dx + real(ws[n])
            y = y₀ + dy + imag(ws[n])
            z = f_base[n].z  # same z position as "base" vortex
            f[n] = (x, y, z)
        end
        update_coefficients!(f)  # needed before any interpolations, plotting, ...
        push!(fs, f)  # add new vortex to list
    end

    # plot_filaments(fs, params.Ls)

    # Mean inter-vortex distance (assuming nearly straight vortices => weak perturbation)
    ℓ = sqrt(params.Ls[1] * params.Ls[2] / length(fs))

    # Rotation rate (Feynman's rule)
    Ω = params.Γ / (2 * ℓ^2)

    #== Run simulation ==#
    # Initialise problem
    # Tsim = 1 * BiotSavart.kelvin_wave_period(params, Lz)  # total simulation time
    Tsim = nrotations * 2π / Ω  # simulate 4 rotation periods
    prob = VortexFilamentProblem(fs, Tsim, params)
    println(prob, '\n')

    # Define callback function to be called at each iteration
    times = Float64[]            # will contain the time associated to each timestep
    energy = Float64[]           # will contain the kinetic energy at each timestep
    vortex_length = Float64[]    # will contain the total vortex length at each timestep
    positions = ComplexF64[]  # will contain vortex positions (w = x + im * y) at each timestpp

    tsf = TimeSeriesFile()

    t_save_step = (prob.tspan[2] - prob.tspan[1]) / (nsave + 1)
    t_save_next = Ref(prob.tspan[1])

    function callback(iter)
        local (; nstep, t,) = iter
        local (; tspan,) = iter.prob
        local (; Ls,) = iter.prob.p
        if nstep == 0  # make sure vectors are empty at the beginning of the simulation
            foreach(empty!, (times, energy, vortex_length, positions, tsf))
            t_save_next[] = tspan[1]
        end
        local E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
        local Lvort = Diagnostics.filament_length(iter; quad = GaussLegendre(3))
        t_rel = (t - tspan[1]) / (tspan[2] - tspan[1])
        if nstep % 10 == 0
            @show nstep, t, t_rel, E, Lvort
        end
        if t ≥ t_save_next[]
            t_save_next[] += t_save_step
            local fname = "rotation_$nstep.vtkhdf"
            save_checkpoint(fname, iter; refinement = 4) do io
                io["curvature"] = CurvatureScalar()
                io["velocity"] = iter.vL
            end
            tsf[t] = fname
            save("rotation.vtkhdf.series", tsf)
        end
        push!(times, t)
        push!(energy, E)
        push!(vortex_length, Lvort)
        # Save locations of the first filament in complex form.
        resize!(positions, length(positions) + N)
        positions_t = @view positions[(end - N + 1):end]  # values to be filled
        interpolate_at_equidistant_zs!(positions_t, iter.fs[1]; Lz = Ls[3])
        nothing
    end

    # Determine timestep and temporal scheme
    δ = Lz / N  # line resolution (assuming nearly straight lines)
    dt_kw = BiotSavart.kelvin_wave_period(params, δ)
    scheme = RK4()
    dt = dt_factor * dt_kw
    # scheme = Strang(RK4(); nsubsteps = 4)  # Strang splitting allows to use larger timesteps
    # dt = 2 * dt_factor * scheme.nsubsteps * dt_kw

    # Initialise and run simulation
    iter = init(prob, scheme; dt, callback)
    println(iter, '\t')
    solve!(iter)
    println(iter.to)  # show timing information

    # Check average evolution of vortex length (wrt length of equilibrium position)
    L_equilibrium = length(fs) * Lz
    vortex_length_extra = @. (vortex_length / L_equilibrium - 1)
    scatterlines(times, vortex_length_extra; axis = (xlabel = "Time", ylabel = "Normalised perturbation amplitude"))
    @show std(vortex_length) / L_equilibrium
    @show (vortex_length[end] - vortex_length[begin]) / L_equilibrium

    (;
        N, Ω, ℓ, params, times, energy, vortex_length, positions,
        w_equilibrium, timer = iter.to,
    )
end

function analyse_kelvin_waves(results)
    (; N, positions, times, params, Ω, ℓ, w_equilibrium) = results

    ws_mat = reshape(positions, N, :)  # reinterpret 1D vector as a 2D matrix: ws_mat[i, j] = w(z_i, t_j)

    t_inds = eachindex(times)[begin:end]  # this may be modified to use a subset of the simulation time
    Nt = length(t_inds)       # number of timesteps to analyse
    ws_h = ws_mat[:, t_inds]  # create copy of positions to avoid modifying the simulation output

    # Remove equilibrium position
    # ws_h .-= w_equilibrium[1]  # 1 = index of vortex

    # Apply window function over temporal dimension
    window = DSP.Windows.hanning(Nt)
    for i in axes(ws_h, 1)
        @. @views ws_h[i, :] *= window
    end

    # Perform a 2D FFT and plot the results
    # fft!(ws_h, (1, 2))
    fft!(ws_h, (1,))  # spatial FFT
    ws_h[1, :] .= 0   # remove spatial mean
    fft!(ws_h, (2,))  # temporal FFT
    ws_h ./= length(ws_h)  # normalise FFT

    # Combine +k and -k wavenumbers and plot as a function of |k| (results should be symmetric anyways)
    ws_abs2 = @views abs2.(ws_h[1:((end + 1) ÷ 2), :])  # from k = 0 to k = +kmax
    @views ws_abs2[2:end, :] .+= abs2.(ws_h[end:-1:(end ÷ 2 + 2), :])  # from -kmax to -1

    Δt = times[t_inds[2]] - times[t_inds[1]]      # timestep
    ks = fftfreq(N, 2π * N / Lz)[1:((end + 1) ÷ 2)]  # wavenumbers (≥ 0 only)
    ωs = fftfreq(Nt, 2π / Δt)     # frequencies
    ωs_shift = fftshift(ωs)
    w_plot = fftshift(log10.(ws_abs2), (2,))  # FFT shift along second dimension (frequencies)

    wplot_max = ceil(Int, maximum(@view w_plot[2:end, :]))

    k_max = ks[end]           # maximum wavenumber
    ω_max = -ωs_shift[begin]  # maximum frequency

    # Analytical dispersion relation according to Raja Gopal (1964) for comparison.
    # This comes from Donnelly 1991, eqs. (6.4)-(6.5).
    ks_fine = range(0, k_max; step = ks[2] / 4)

    # Estimate J(kℓ) = ∫_{0}^{kℓ} K₀(s) s ds = ℓ² ∫_{0}^{k} K₀(qℓ) q dq for each k ∈ ks_fine.
    Js = let
        Js = similar(ks_fine)
        @assert ks_fine[begin] == 0
        integrate_incremental!(Js, ks_fine; quad = GaussLegendre(4)) do q
            qℓ = q * ℓ
            iszero(qℓ) ? qℓ : ℓ^2 * besselk0(qℓ) * q  # avoid fake divergence at qℓ = 0 (the function is actually well defined)
        end
        Js
    end

    ωs_rajagopal = let
        local (; Γ, a, Δ) = params
        local γ = MathConstants.eulergamma
        @. -(2 - Js) * Ω - Γ * ks_fine^2 / (4 * π) * (
            # log(2 / (abs(ks_fine) * a)) - γ + 1/2 - Δ
            log(2 / (abs(ks_fine) * a)) + 1 / 2 - Δ
        )
    end

    # Small kℓ limit -> J = 0
    ωs_small_k = let
        local (; Γ, a, Δ) = params
        local γ = MathConstants.eulergamma
        @. -2Ω - Γ * ks_fine^2 / (4 * π) * (
            # log(2 / (abs(ks_fine) * a)) - γ + 1/2 - Δ
            log(2π / (abs(ks_fine) * a)) + 1 / 2 - Δ
        )
    end

    # Small kℓ: 2π/k is replaced by ℓ in the logarithm.
    ωs_ℓ = let
        local (; Γ, a, Δ) = params
        @. -2Ω - Γ * ks_fine^2 / (4 * π) * (log(ℓ / a) + 1 / 2 - Δ)
    end

    k_ℓ = 2π / ℓ  # wavenumber associated to inter-vortex distance

    fontsize = 14
    fig = Figure(; size = (400, 300), fontsize)
    ax = Axis(fig[1, 1]; xlabel = L"kℓ", ylabel = L"ω(k)/Ω")
    ax.xgridvisible = ax.ygridvisible = false
    ax.xlabelpadding = 4
    ax.ylabelpadding = 0
    xlims!(ax, 0.9 * k_max .* (0, 1))
    # ylims!(ax, 0.6 * ω_max .* (-1, 0))
    ylims!(ax, (-4, 0))
    hm = @views heatmap!(
        ax, ks[2:end] .* ℓ, ωs_shift ./ Ω, w_plot[2:end, :];
        rasterize = 8,
        # colormap = Reverse(:deep),
        # colormap = :OrRd,
        colormap = Reverse(:oslo),
        # colormap = Reverse(:afmhot),
        # colorrange = round.(Int, extrema(w_plot[2:end, :])) .+ (3, 0),
        colorrange = (-16, -2),
        # colorrange = (wplot_max - 3, wplot_max - 2),
        # colorrange = (wplot_max - 3, wplot_max - 2),
    )
    Colorbar(fig[1, 2], hm; label = L"\log \, |\hat{w}|^2")
    hlines!(ax, 0; color = :black, linestyle = :dash, linewidth = 1)
    let color = :chartreuse3
        ls = hlines!(ax, -2; color, linestyle = :dash, linewidth = 2)
        pos_end = @lift($(ls.points)[2])  # the hline connects two points: [start, end] (in data coordinates)
        text!(ax, pos_end; text = L"-2Ω", offset = (-4, -4), align = (:right, :top), color, fontsize = fontsize + 2)
    end
    # let color = :lightblue
    #     vlines!(ax, k_ℓ; color, linestyle = :dot)
    #     text!(ax, @lift((k_ℓ, $ybottom)); text = L"\frac{2π}{ℓ}", offset = (6, 4), align = (:left, :bottom), color)
    # end
    # lines!(ax, ks_fine, ωs_small_k; color = :orangered, linestyle = :dash, linewidth = 1.5)
    # lines!(ax, ks_fine, ωs_rajagopal; color = (:lightgreen, 0.8), linestyle = :dash, linewidth = 1.5, label = L"(2 - J) Ω + \frac{κ}{4π} \ln(2π/ka)")
    let color = :chocolate, n = searchsortedlast(ks_fine, 6.0)
        lines!(ax, ks_fine .* ℓ, ωs_ℓ ./ Ω; color, linestyle = :dash, linewidth = 2.0, label = L"2Ω + \frac{κ k^2}{4π} \ln(ℓ/a)")
        annotation!(
            ax,
            2 * ks_fine[n] * ℓ, -3.5,     # label position
            ks_fine[n] * ℓ, ωs_ℓ[n] / Ω;  # target position
            labelspace = :data,
            color, text = L"-2Ω - \frac{κ k^2}{4π} \ln(ℓ/a)", fontsize = fontsize + 2,
            style = Ann.Styles.LineArrow(),
        )
    end
    # axislegend(ax; position = (0, 0), framevisible = false, labelcolor = :white, labelsize = 16)
    # DataInspector(fig)
    xlims!(ax, 0, 1.25)
    ax.xticks = range(0, 1.2; step = 0.2)
    # ylims!(ax, -120, 0)
    let bg = fig.scene.backgroundcolor[]
        fig.scene.backgroundcolor[] = RGBAf(0, 0, 0, 0)
        # save("rotation_dispersion.svg", fig; backend = CairoMakie)
        fig.scene.backgroundcolor[] = bg
    end
    fig
end

## Run simulation

Lz = 2π / 1   # domain period (vertical)
Lx = 2π / 32  # domain period (horizontal) -- this is proportional to the inter-vortex distance

Ly = Lx * sqrt(3) / 2  # hexagonal lattice (this reduces aspect ratio Lx:Ly)
nvortices = (2 * 1)^2  # must be (2n)² for aspect ratio 2:√3

Γ = 1.0   # circulation
a = 1e-7  # vortex core size
ℓ = sqrt(Lx * Ly / nvortices)
Ω = Γ / (2 * ℓ^2)
kz_min = 2π / Lz
(Γ / 4π) * kz_min^2 * log(ℓ / a)

# Initial vortices
N = 64            # number of discretisation points per line
δ = Lz / N
@show δ / ℓ

in_phase_perturbation = true  # perturbations of all vortices are in phase?
Az_rms = 1e-4 * Lz   # amplitude of random perturbation (rms value)
Ax = 0e-2 * Lx
Ay = Ax
Axy_rms = 0 * ℓ / 100  # amplitude of random horizontal perturbation
kx = 1; ky = 1
nrotations = 4
nsave = 100 * nrotations
rcut_factor = 0.4  # for performance: balance between short and long-range (max 2/5)
dt_factor = 0.5    # small dt probably needed for strong rotations

# Test imaginary time evolution
# let Ly = Lx * sqrt(3), rng = StableRNG(55), Nvort = 64
#     local xs = Lx * rand(rng, Nvort)
#     local ys = Ly * rand(rng, Nvort)
#     find_stable_equilibrium!(xs, ys; Ls = (Lx, Ly, Lz), verbose = true, checkpoints = true)
# end

method = QuinticSplineMethod()  # filament discretisation method

results = run_simulation(;
    nvortices, Γ, a, Ls = (Lx, Ly, Lz), N, in_phase_perturbation,
    Az_rms, method,
    Ax, Ay, kx, ky, nrotations,
    Axy_rms,
    rcut_factor, dt_factor,
)

## Spatiotemporal analysis

analyse_kelvin_waves(results)
