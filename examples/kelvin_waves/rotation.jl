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
using GLMakie  # for plots
using InverseFunctions: InverseFunctions

## Helper functions

function init_biot_savart_parameters(;
        Ls, β = 3.5,
        rcut = min(Ls...) * (2/5),  # default = maximum possible cut-off distance for short-range part (depends on CellListsBackend parameters)
        kws...,
    )
    α = β / rcut                # Ewald splitting parameter
    kmax = 2 * α * β            # maximum resolved wavenumber (Nyquist frequency) for long-range part
    Ns = map(Ls) do L
        1 + ceil(Int, kmax * L / π)  # long-range grid resolution (~FFT size)
    end
    ParamsBiotSavart(;
        Γ = 1.0,    # vortex circulation
        a = 1e-8,   # vortex core size
        Δ = 1/4,    # vortex core parameter (1/4 for a constant vorticity distribution)
        α = α,      # Ewald splitting parameter
        Ls,  # same domain size in all directions
        Ns,  # same long-range resolution in all directions
        rcut = β / α,    # cut-off distance for short-range computations
        quadrature = GaussLegendre(3),        # quadrature for integrals over filament segments
        quadrature_near_singularity = GaussLegendre(3),
        lia_segment_fraction = 0.2,
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4)),  # this is the default
        backend_short = CellListsBackend(2),
    )
end

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
function random_perturbation!(rng, ws::AbstractVector{<:Complex}; Lz::Real, A_rms::Real)
    N = length(ws)

    # First define perturbation in Fourier space
    ks = fftfreq(N, 2π * N / Lz)  # wavenumbers
    kmax = π * N / Lz
    kmax_perturb = 0.5 * kmax   # perturb up to this wavenumber (to avoid issues at the discretisation scale)
    for i in eachindex(ws, ks)
        kabs = abs(ks[i])
        if 0 < kabs <= kmax_perturb
            ws[i] = randn(rng, ComplexF64)  # perturb all wavenumbers equally
        else
            ws[i] = 0
        end
    end

    # Rescale perturbation to obtain the wanted amplitude (rms value)
    factor = A_rms / sqrt(sum(abs2, ws))
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

# Adapted from MakieExtra.jl: transform between :data and :relative spaces in Makie plots.
data2rel(ax::Axis, which::Symbol, orig) = transform_val_space(ax, which, :data => :relative, orig)
rel2data(ax::Axis, which::Symbol, orig) = transform_val_space(ax, which, :relative => :data, orig)

transform_val_space(ax, which::Symbol, args...) = transform_val_space(ax, Dict(:x=>1, :y=>2)[which], args...)

function transform_val_space(ax::Axis, which::Int, spaces::Pair{Symbol,Symbol}, orig::Union{Tuple,AbstractVector})
    scene = Makie.get_scene(ax)
    lift(scene.camera.projectionview, Makie.plots(ax)[1].model, Makie.transform_func(ax), scene.viewport, orig) do _, _, tf, _, orig
        tf_cur = Base.setindex(tf, identity, 3 - which)
        if spaces[1] == :data
            # orig = map(|>, orig, tf_cur)
            orig = map((f, x) -> f(x), tf_cur, orig)
        end
        new = Makie.project(scene.camera, spaces..., Point2(orig))
        if spaces[2] == :data
            new = Point2(map((f, x) -> InverseFunctions.inverse(f)(x), tf_cur, new))
        end
        Base.setindex(orig, new[which], which)
    end
end

function transform_val_space(ax::Axis, which::Int, spaces::Pair{Symbol,Symbol}, orig::Number)
    new = transform_val_space(ax, which, spaces, Base.setindex((1, 1), orig, which))
    @lift $new[which]
end

function run_simulation(;
        lattice::Symbol = :square,   # :square or :hexagonal
        Lz = 2π,             # domain period (vertical)
        Lx = Lz,             # domain period (horizontal) -- proportional to inter-vortex distance
        N = 64,              # number of points per vortex (= line resolution)
        in_phase_perturbation = true,  # perturbations of all vortices are in phase?
        nvort_per_dir = 2,   # number of vortices along X direction (number of vortices along Y depends on this and on `lattice` value)
        A_rms = 1e-5 * Lz,   # amplitude of random perturbation (rms value)
        method = QuinticSplineMethod(),  # filament discretisation method
        rcut_factor = 2/5,
    )
    #== Biot-Savart parameters ==#
    # Domain size in y direction
    Ly = if lattice == :square
        Lx
    elseif lattice == :hexagonal
        Lx * sqrt(3)
    else
        throw(ArgumentError("`lattice` parameter should be :square or :hexagonal"))
    end
    Ls = (Lx, Ly, Lz)
    rcut = min(Ls...) * rcut_factor
    params = init_biot_savart_parameters(; Ls, rcut)
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
    w_equilibrium = ComplexF64[]  # equilibrium position (x + iy) of each filament
    ws = zeros(ComplexF64, N)  # this will be used to store random perturbations
    rng = StableRNG(42)        # initialise random number generator (RNG)

    # Relative position of one or more vortices in cell (same for x and y coordinates)
    cell_positions = if lattice == :square
        [1/2]       # 1 vortex in the middle of the cell
    elseif lattice == :hexagonal
        [1/4, 3/4]  # 2 vortices in a minimal periodic unit cell
    end

    for j in 1:nvort_per_dir, i in 1:nvort_per_dir, pos in cell_positions
        f = similar(f_base)
        if !in_phase_perturbation || isempty(fs)
            # Don't update ws if in_phase_perturbation == true and this is not the first vortex.
            random_perturbation!(rng, ws; Lz, A_rms)  # compute random perturbation, modifying ws
        end
        x₀ = (i - 1 + pos) * params.Ls[1] / nvort_per_dir  # x position of vortex (before perturbation)
        y₀ = (j - 1 + pos) * params.Ls[2] / nvort_per_dir  # y position of vortex (before perturbation)
        for n in eachindex(f, f_base)
            x = x₀ + real(ws[n])
            y = y₀ + imag(ws[n])
            z = f_base[n].z  # same z position as "base" vortex
            f[n] = (x, y, z)
        end
        update_coefficients!(f)  # needed before any interpolations, plotting, ...
        push!(fs, f)  # add new vortex to list
        push!(w_equilibrium, x₀ + im * y₀)
    end

    # plot_filaments(fs, params.Ls)

    # Mean inter-vortex distance (assuming nearly straight vortices => weak perturbation)
    ℓ = sqrt(params.Ls[1] * params.Ls[2] / length(fs))  # = Lh / nvort_per_dir

    # Rotation rate (Feynman's rule)
    Ω = params.Γ / (2 * ℓ^2)

    #== Run simulation ==#
    # Initialise problem
    # Tsim = 1 * BiotSavart.kelvin_wave_period(params, Lz)  # total simulation time
    Tsim = 4 * 2π / Ω  # simulate 4 rotation periods
    prob = VortexFilamentProblem(fs, Tsim, params)
    println(prob, '\n')

    # Define callback function to be called at each iteration
    times = Float64[]            # will contain the time associated to each timestep
    energy = Float64[]           # will contain the kinetic energy at each timestep
    vortex_length = Float64[]    # will contain the total vortex length at each timestep
    positions = ComplexF64[]  # will contain vortex positions (w = x + im * y) at each timestpp

    tsf = TimeSeriesFile()

    t_save_step = (prob.tspan[2] - prob.tspan[1]) / 101
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
    dt = dt_kw
    # scheme = Strang(RK4(); nsubsteps = 4)  # Strang splitting allows to use larger timesteps
    # dt = 2 * scheme.nsubsteps * dt_kw

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
        w_equilibrium,
    )
end

## Run simulation

lattice = :square
# lattice = :hexagonal
Lz = 2π       # domain period (vertical)
Lx = Lz / 8   # domain period (horizontal) -- this is proportional to the inter-vortex distance

# Initial vortices
N = 64             # number of discretisation points per line
in_phase_perturbation = false  # perturbations of all vortices are in phase?
nvort_per_dir = 4   # number of vortices per direction (x, y)
A_rms = 1e-5 * Lz   # amplitude of random perturbation (rms value)
method = QuinticSplineMethod()  # filament discretisation method

results = run_simulation(;
    lattice, Lz, Lx, N, in_phase_perturbation,
    nvort_per_dir, A_rms, method,
)

## Spatiotemporal analysis

(; N, positions, times, params, Ω, ℓ, w_equilibrium,) = results

ws_mat = reshape(positions, N, :)  # reinterpret 1D vector as a 2D matrix: ws_mat[i, j] = w(z_i, t_j)

t_inds = eachindex(times)[begin:end]  # this may be modified to use a subset of the simulation time
Nt = length(t_inds)       # number of timesteps to analyse
ws_h = ws_mat[:, t_inds]  # create copy of positions to avoid modifying the simulation output

# Remove equilibrium position
ws_h .-= w_equilibrium[1]  # 1 = index of vortex

# Apply window function over temporal dimension
window = DSP.Windows.hanning(Nt)
for i ∈ axes(ws_h, 1)
    @. @views ws_h[i, :] *= window
end

# Perform a 2D FFT and plot the results
fft!(ws_h, (1, 2))
ws_h ./= length(ws_h)  # normalise FFT

# Combine +k and -k wavenumbers and plot as a function of |k| (results should be symmetric anyways)
ws_abs2 = @views abs2.(ws_h[1:((end + 1)÷2), :])  # from k = 0 to k = +kmax
@views ws_abs2[2:end, :] .+= abs2.(ws_h[end:-1:(end÷2 + 2), :])  # from -kmax to -1

Δt = times[t_inds[2]] - times[t_inds[1]]      # timestep
ks = fftfreq(N, 2π * N / Lz)[1:(end + 1)÷2]  # wavenumbers (≥ 0 only)
ωs = fftfreq(Nt, 2π / Δt)     # frequencies
ωs_shift = fftshift(ωs)
w_plot = fftshift(log10.(ws_abs2), (2,))  # FFT shift along second dimension (frequencies)

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
    local (; Γ, a, Δ,) = params
    local γ = MathConstants.eulergamma
    @. -(2 - Js) * Ω - Γ * ks_fine^2 / (4 * π) * (
        # log(2 / (abs(ks_fine) * a)) - γ + 1/2 - Δ
        log(2π / (abs(ks_fine) * a)) + 1/2 - Δ
    )
end

# Small kℓ limit -> J = 0
ωs_small_k = let
    local (; Γ, a, Δ,) = params
    local γ = MathConstants.eulergamma
    @. -2Ω - Γ * ks_fine^2 / (4 * π) * (
        # log(2 / (abs(ks_fine) * a)) - γ + 1/2 - Δ
        log(2π / (abs(ks_fine) * a)) + 1/2 - Δ
    )
end

# Small kℓ: 2π/k is replaced by ℓ in the logarithm.
ωs_ℓ = let
    local (; Γ, a, Δ,) = params
    @. -2Ω - Γ * ks_fine^2 / (4 * π) * (log(ℓ / a) + 1/2 - Δ)
end

k_ℓ = 2π / ℓ  # wavenumber associated to inter-vortex distance

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"k", ylabel = L"ω", xlabelsize = 20, ylabelsize = 20)
xlims!(ax, 0.8 * k_max .* (0, 1))
ylims!(ax, 0.8 * ω_max .* (-1, 1))
hm = heatmap!(
    ax, ks, ωs_shift, w_plot;
    colormap = Reverse(:deep),
    colorrange = round.(Int, extrema(w_plot)) .+ (3, 0),
)
Colorbar(fig[1, 2], hm; label = L"\log \, |\hat{w}|^2", labelsize = 20)
hlines!(ax, 0; color = :white, linestyle = :dash, linewidth = 1)
xright = rel2data(ax, :x, 1.0)  # rightmost coordinate (in data space)
ybottom = rel2data(ax, :y, 0.0)
text!(ax, L"0"; position = @lift(($xright, 0)), offset = (-4, +2), align = (:right, :bottom), color = :white, fontsize = 16)
hlines!(ax, -2 * Ω; color = :yellow, linestyle = :dash, linewidth = 1)
text!(ax, L"-2Ω"; position = @lift(($xright, -2 * Ω)), offset = (-4, -2), align = (:right, :top), color = :yellow, fontsize = 16)
let color = :lightblue
    vlines!(ax, k_ℓ; color, linestyle = :dot)
    text!(ax, L"\frac{2π}{ℓ}"; position = @lift((k_ℓ, $ybottom)), offset = (6, 4), align = (:left, :bottom), color, fontsize = 16)
end
# lines!(ax, ks_fine, ωs_small_k; color = :orangered, linestyle = :dash, linewidth = 1.5)
lines!(ax, ks_fine, ωs_rajagopal; color = (:lightgreen, 0.8), linestyle = :dash, linewidth = 1.5, label = L"(2 - J) Ω + \frac{κ}{4π} \ln(2π/ka)")
lines!(ax, ks_fine, ωs_ℓ; color = (:yellow, 0.8), linestyle = :dash, linewidth = 2.0, label = L"2Ω + \frac{κ}{4π} \ln(ℓ/a)")
# axislegend(ax; position = (0, 0), framevisible = false, labelcolor = :white, labelsize = 16)
# DataInspector(fig)
fig
