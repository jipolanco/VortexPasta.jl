# Simple vortex ring example in a periodic box.
# This script can be used to test the stability of different timestepping methods.

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

function generate_ring_filament(
        ::Type{T}, Np, method = QuinticSplineMethod();
        R, Ls,
    ) where {T}
    S = define_curve(Ring(); scale = R, translate = Ls ./ 2)
    [Filaments.init(S, ClosedFilament{T}, Np, method)]
end

function run_vortex_ring(
        params::ParamsBiotSavart{T};
        method = QuinticSplineMethod(),
        dissipation = NoDissipation(),
        dt_factor = 0.5,
        scheme = RK4(),
        outdir = "output",
        Np = 64,
        R = T(1.0),
        tsim_factor = 1,
    ) where {T}
    # Initialise problem
    fs = generate_ring_filament(T, Np, method; R, Ls = params.Ls)
    tsim = tsim_factor * R^2 / params.Γ
    prob = VortexFilamentProblem(fs, tsim, params)
    println(prob)

    # Initialise solver
    δ_min = minimum_node_distance(fs)
    dt = dt_factor * BiotSavart.kelvin_wave_period(params, δ_min)

    mkpath(outdir)

    times = T[]
    vortex_length = T[]
    energies = T[]
    impulse = Vec3{T}[]
    tsf = TimeSeriesFile()

    function callback(iter::VortexFilamentSolver)
        local (; nstep, t, dt, fs, to,) = iter
        local quad = params.quad
        local E = Diagnostics.kinetic_energy(iter; quad)
        local Lvort = Diagnostics.filament_length(iter; quad)
        local p⃗ = Diagnostics.vortex_impulse(iter; quad)

        push!(times, t)
        push!(energies, E)
        push!(vortex_length, Lvort)
        push!(impulse, p⃗)

        if nstep % 10 == 0
            cd(outdir) do
                local fname = "vortex_ring_$nstep.vtkhdf"
                write_vtkhdf(fname, fs; refinement = 1) do io
                    io["velocity"] = iter.vs
                    io["Curvature"] = CurvatureScalar()
                end
                tsf[t] = fname
                save("vortex_ring.vtkhdf.series", tsf)
            end
        end

        if nstep % 10 == 0
            let E₀ = energies[begin]
                @show nstep, t/tsim, dt, E/E₀
            end
        end
    end

    iter = init(
        prob, scheme;
        dt, dissipation, fold_periodic = false,
        callback,
    )
    println(iter)
    solve!(iter)

    println(iter.to)

    (; iter, params, prob, times, energies, impulse, vortex_length,)
end

##

Ls = (2π, 2π, 2π)
Ns = (32, 32, 32)
params = generate_biot_savart_parameters(; Ls, Ns)

dissipation = NoDissipation()
# dissipation = SmallScaleDissipationBS(; kdiss = 14, α = 1.0)

outdir = "output"

results = run_vortex_ring(
    params;
    method = QuinticSplineMethod(),
    dt_factor = 1.0,
    dissipation, outdir,
    scheme = SSPRK33(),
    tsim_factor = 1.0,
    Np = 64, R = 1.0,
);

##

(; iter, times, energies, impulse, vortex_length,) = results

impulse_norm = map(norm, impulse);

fig = Figure()
ax = Axis(fig[1, 1]; xlabel = "Time", ylabel = "Relative change")
# ylims!(ax, 0.995, 1.001)
scatterlines!(ax, times, energies ./ first(energies); label = "Energy")
scatterlines!(ax, times, impulse_norm ./ first(impulse_norm); label = "Impulse")
scatterlines!(ax, times, vortex_length ./ first(vortex_length); label = "Vortex length")
axislegend(ax; position = (0, 0))
fig
