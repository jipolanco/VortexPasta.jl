# This is adapted from the background_vorticity.jl (rotation frame) test.

using VortexPasta
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.PredefinedCurves
using VortexPasta.FilamentIO
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using Test

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

function init_params_biot_savart(; L, β)
    rcut = L * 2 / 5
    α = β / rcut
    kmax = 2α * β
    M = floor(Int, L * kmax / π)
    ParamsBiotSavart(
        Float64;
        Γ = 1.0, a = 1e-8,
        α, Ls = L,
        rcut, Ns = (M, M, M),
        backend_short = CellListsBackend(2),
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4))  # ~1e-6 accuracy
    )
end

function test_with_rotation(; save_outputs = false,)
    L = 2π
    β = 3.5
    params = init_params_biot_savart(; L, β)

    Np = 16
    method = QuinticSplineMethod()

    # Create 3×3 array of aligned vortices with identical helical perturbations.
    # The 0.98 factor means that vortices are not exactly homogeneously distributed at the
    # beginning, and they should move until they asymptotically form a regular lattice.
    Nx = 3
    xs = range(0, 0.98 * L; length = 2 * Nx + 1)[2:2:end]
    ys = xs

    fs = map(Iterators.product(xs, ys)) do (x, y)
        p = PeriodicLine(r = t -> cispi(4t) / 100)
        S = define_curve(p; scale = L, translate = (x, y, L / 2), orientation = +1)
        ζs = range(0, 1; length = 2 * Np + 1)[2:2:end]
        Filaments.init(S, ClosedFilament, ζs, method)
    end |> vec

    Tsim = 10
    prob = VortexFilamentProblem(fs, Tsim, params)

    vortex_length_equilibrium = length(fs) * L  # total vortex length at equilibrium
    vortex_distance_equilibrium = L / Nx  # minimal distance between two vortices at equilibrium

    iter = init(
        prob, RK4();
        dt = 1.0, adaptivity = AdaptBasedOnSegmentLength(0.5),
        mode = MinimalEnergy(),
    )

    vortex_length_err = Inf    # this is 0 if all vortices are perfectly straight (converges fast)
    vortex_distance_err = Inf  # this is 0 if vortices form a perfectly regular lattice (takes more time to converge)

    energies = Float64[]
    tsf = TimeSeriesFile()

    while iter.t < 10.0
        if save_outputs && iter.nstep % 10 == 0
            filename = "min_energy_$(iter.nstep).vtkhdf"
            save_checkpoint(filename, iter) do io
                io["velocity"] = iter.vs
            end
            tsf[iter.t] = filename
            save("min_energy.vtkhdf.series", tsf)
        end
        step!(iter)
        (; nstep, t) = iter
        Lvort = Diagnostics.filament_length(iter; quad = GaussLegendre(3))
        E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
        v²_max = maximum(iter.vs) do vs
            maximum(v -> sum(abs2, v), vs)
        end
        vortex_distance = (iter.fs[2][1] - iter.fs[1][1]).x
        vortex_length_err = abs(Lvort - vortex_length_equilibrium) / vortex_length_equilibrium
        vortex_distance_err = abs(vortex_distance - vortex_distance_equilibrium) / vortex_distance_equilibrium
        push!(energies, E)
        if VERBOSE && iter.nstep % 10 == 0
            @show nstep, t, vortex_length_err, E, v²_max / E, vortex_distance_err
        end
        vortex_distance_err < 0.01 && break
    end

    @test all(<(0), diff(energies))   # energy always decreases in time
    @test vortex_length_err == 0      # perturbations were completely dissipated (vortices are straight)
    @test vortex_distance_err < 0.01  # vortices tend to form a regular lattice (slowly...)

    nothing
end

@testset "Minimal energy mode" begin
    @testset "Rotation" test_with_rotation()
end
