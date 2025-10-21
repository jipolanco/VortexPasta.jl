# Reconnection of two vortex rings initially on the same plane.

using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.BiotSavart
using VortexPasta.FilamentIO
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using VortexPasta.Forcing
using VortexPasta
using LinearAlgebra: norm, ⋅
using Rotations: Rotations
using TimerOutputs: @timeit
using GLMakie

function generate_biot_savart_parameters(;
        Ls::NTuple{3, T}, β = 3.5,
        Ns::Union{Nothing, Dims{3}} = nothing,
        rcut::Union{Nothing, Real} = Ns === nothing ? min(Ls...) / 3 : nothing,
    ) where {T}
    Γ = 1.0
    a = 1e-7
    Δ = 1/2
    if Ns === nothing
        @assert rcut isa Real
        α = β / rcut
        kmax = 2α * β
        Ns = ceil.(Int, (kmax / π) .* Ls)
    else
        @assert rcut === nothing
        ks_max = @. π * (Ns - 1) / Ls
        kmax = min(ks_max...)
        α = kmax / (2 * β)
        rcut = β / α
    end
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = NonuniformFFTsBackend(CPU(); σ = T(1.5), m = HalfSupport(4)),
        # backend_long = NonuniformFFTsBackend(ROCBackend(); σ = T(1.5), m = HalfSupport(4), gpu_method = :shared_memory),
        quadrature = GaussLegendre(3),
        quadrature_near_singularity = GaussLegendre(3),
        lia_segment_fraction = 0.2,
    )
end

function generate_ring_filaments(
        ::Type{T}, Np, method = QuinticSplineMethod();
        R, distance,
        angle_degrees = 0.0,
    ) where {T}
    p = Ring()
    h = distance / 2
    orientation = -1
    scale = (R, R, R)
    α = deg2rad(angle_degrees) / 2
    rot_l = Rotations.RotX(-α)  # rotation of -α about the X axis
    rot_r = Rotations.RotX(+α)  # rotation of +α about the X axis
    Sl = define_curve(p; scale, orientation, rotate = rot_l, translate = (-(R + h), 0, 0))  # left ring
    Sr = define_curve(p; scale, orientation, rotate = rot_r, translate = (+(R + h), 0, 0))  # right ring
    [
        Filaments.init(Sl, ClosedFilament{T}, Np, method),
        Filaments.init(Sr, ClosedFilament{T}, Np, method),
    ]
end

function run_reconnections(
        params::ParamsBiotSavart{T};
        method = QuinticSplineMethod(),
        dissipation = NoDissipation(),
        dt_factor = 0.5,
        scheme = RK4(),
        outdir = "output",
        Np = 256,
        R = T(1.0),
        distance = R / 40,
        tsim_factor = 1,
    ) where {T}
    # Initialise vortex rings
    angle_degrees = 0

    fs = generate_ring_filaments(T, Np, method; R, distance, angle_degrees,)

    # Initialise problem
    tsim = tsim_factor * distance^2 / params.Γ
    prob = VortexFilamentProblem(fs, tsim, params)
    println(prob)

    # Initialise solver
    δ = minimum_node_distance(fs)
    δ_min = 0.8 * δ
    reconnect = ReconnectFast(δ_min)
    refinement = RefineBasedOnSegmentLength(δ_min)
    dt = dt_factor * BiotSavart.kelvin_wave_period(params, δ_min)

    # adaptivity = AdaptBasedOnSegmentLength(dt_factor)
    adaptivity = NoAdaptivity()

    iter = init(
        prob, scheme;
        dt,
        reconnect,
        dissipation,
        refinement, adaptivity,
        fold_periodic = false,
    )

    mkpath(outdir)

    # Run solver
    times = T[]
    dts = T[]
    vortex_length = T[]
    energies = T[]
    impulse = Vec3{T}[]
    helicities = T[]

    tsf = TimeSeriesFile()

    while iter.t < tsim
        local (; nstep, t, dt, fs, to,) = iter
        local quad = params.quad
        @timeit to "Energy" E = Diagnostics.kinetic_energy(iter; quad)
        Lvort = Diagnostics.filament_length(iter; quad)
        p⃗ = Diagnostics.vortex_impulse(iter; quad)
        @timeit to "Helicity" H = Diagnostics.helicity(iter; quad)
        push!(times, t)
        push!(energies, E)
        push!(helicities, H)
        push!(vortex_length, Lvort)
        push!(impulse, p⃗)
        push!(dts, dt)

        cd(outdir) do
            local fname = "reconnect_points_$nstep.vtkhdf"
            @timeit to "Write VTKHDF (R1)" write_vtkhdf(fname, fs; refinement = 1) do io
                io["velocity"] = iter.vs
                io["Curvature"] = CurvatureScalar()
            end
            tsf[t] = fname
            save("reconnect_points.vtkhdf.series", tsf)
        end

        if nstep % 10 == 0
            @show nstep, t/tsim, dt, length(fs), E
        end
        # length(fs) == 1 && break  # stop right after reconnection
        step!(iter)
    end

    println(iter.to)

    (; iter, params, prob, times, dts, energies, helicities, impulse, vortex_length,)
end

##

Ls = (2π, 2π, 2π)
Ns = (32, 32, 32)
params = generate_biot_savart_parameters(; Ls, Ns)

dissipation = NoDissipation()
# dissipation = SmallScaleDissipationBS(; kdiss = 14, α = 1.0)

outdir = "output"

results = run_reconnections(
    params;
    method = QuinticSplineMethod(),
    dt_factor = 1.0,
    dissipation, outdir,
    scheme = RK4(),
    tsim_factor = 10.0,
);

##

(; iter, times, energies, helicities, impulse, vortex_length,) = results

# quad = GaussLegendre(4)
# @btime Diagnostics.kinetic_energy($iter; quad = $quad)

impulse_norm = map(norm, impulse);

fig = Figure()
ax = Axis(
    fig[1, 1]; xlabel = "Time", ylabel = "Relative change",
    # xticklabelsize = 20,
    # yticklabelsize = 20,
)
ylims!(ax, 0.995, 1.001)
scatterlines!(ax, times, energies ./ first(energies); label = "Energy")
scatterlines!(ax, times, impulse_norm ./ first(impulse_norm); label = "Impulse")
# scatterlines!(ax, times, vortex_length ./ first(vortex_length); label = "Vortex length")
axislegend(ax; position = (0, 0))
fig

save("energy.png", fig)
