# Adapted from benchmark/benchmarks.jl

using VortexPasta
using VortexPasta.PredefinedCurves
using VortexPasta.VectorsOfVectors: VectorOfVectors
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Reconnections
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using VortexPasta.Forcing
using Rotations: Rotations
using StableRNGs: StableRNG
using BenchmarkTools
using TimerOutputs
using HDF5: HDF5

using ThreadPinning

if haskey(ENV, "SLURM_JOB_ID")
    pinthreads(:affinitymask)
else
    pinthreads(:cores)
end
threadinfo()

function generate_biot_savart_parameters(
        splitting;
        backend_short = CellListsBackend(CPU(), 2),
        backend_long = NonuniformFFTsBackend(CPU(); σ = 1.5, m = HalfSupport(4)),
        kws...,
    )
    ParamsBiotSavart(;
        Γ = 1.0, a = 1.0e-7, Δ = 1 / 2, splitting,
        backend_short, backend_long,
        quadrature = GaussLegendre(3),
        quadrature_near_singularity = GaussLegendre(3),
        lia_segment_fraction = 0.2,
        kws...,
    )
end

abstract type VFMInitialCondition end

@kwdef struct RandomRings{T} <: VFMInitialCondition
    nrings::Int
    R_min::T
    R_max::T
end

function generate_filaments(rings::RandomRings; Ls, l_res, method)
    (; nrings, R_min, R_max) = rings
    rng = StableRNG(42)
    map(1:nrings) do _
        T = eltype(Ls)
        R = R_min + rand(rng) * (R_max - R_min)
        N = ceil(Int, 2π * R / l_res)
        translate = Ls .* rand(rng, Vec3{T})
        rotate = rand(rng, Rotations.QuatRotation)
        S = define_curve(Ring(); scale = R, translate, rotate)
        Filaments.init(S, ClosedFilament{T}, N, method)
    end
end

function run_benchmark(Np_wanted; ρ_grid = 0.1, Np_per_ring = 1000, params_kws...)
    nrings = Np_wanted ÷ Np_per_ring
    if nrings < 1
        error("Np_wanted should be larger than Np_per_ring")
    end
    Np = nrings * Np_per_ring

    L = 2π
    Ls = (L, L, L)
    R_ring = 0.5
    l_res = 2π * R_ring / Np_per_ring  # typical line resolution
    method = QuinticSplineMethod()
    initial_condition = RandomRings(nrings = nrings, R_min = R_ring, R_max = R_ring)  # generate circular rings
    fs = VectorOfVectors(generate_filaments(initial_condition; Ls, l_res, method))
    @assert sum(length, fs) == Np

    Ngrid = round(Int, cbrt(Np / ρ_grid))
    Ns = (Ngrid, Ngrid, Ngrid)
    splitting = KaiserBesselSplitting(; Ls, Ns, β = 14.0)
    params = generate_biot_savart_parameters(splitting; params_kws...)
    cache = BiotSavart.init_cache(params)

    @info "Running benchmark" Np params

    # fields = (; velocity = similar(fs), streamfunction = similar(fs))  # compute velocity and streamfunction
    fields = (; velocity = similar(fs),)  # compute velocity only
    BiotSavart.compute_on_nodes!(fields, cache, fs)  # warmup (probably not needed with BenchmarkTools...)
    bench = @benchmarkable BiotSavart.compute_on_nodes!($fields, $cache, $fs) evals=1 samples=10 gcsample=true seconds=3600  # 10 runs (or less if we reach 1 hour)
    reset_timer!(cache.to)
    results = run(bench)::BenchmarkTools.Trial

    t_shortrange = let key = "Short-range component (async)"
        to = cache.to[key]
        TimerOutputs.time(to) / TimerOutputs.ncalls(to) / 1e9
    end
    t_longrange = let key = "Long-range component (async)"
        to = cache.to[key]
        TimerOutputs.time(to) / TimerOutputs.ncalls(to) / 1e9
    end

    println(cache.to)

    (; bench = results, timer = cache.to, Np, ρ_grid, params, t_shortrange, t_longrange)
end

function run_benchmark_series(
        Np_wanted_all::AbstractVector;
        ρ_grid = 0.1, outdir = "results",
        backend_short = CellListsBackend(CPU(), 2),
        backend_long = NonuniformFFTsBackend(CPU(); σ = 1.5, m = HalfSupport(4)),
    )
    @info "Running benchmark series" Np_wanted_all ρ_grid backend_short backend_long
    results_all = map(Np_wanted_all) do Np_float
        Np_wanted = round(Int, Np_float)
        println(stdout)
        println(stderr)
        GC.gc()
        results = @time run_benchmark(Np_wanted; ρ_grid, backend_short, backend_long)
        show(stdout, MIME"text/plain"(), results.bench)
        flush(stdout)
        flush(stderr)
        results
    end
    Np_actual = getproperty.(results_all, :Np)
    times_min = map(results_all) do res
        minimum(res.bench.times) / 1e9
    end
    outfile = joinpath(outdir, "results_rho$(ρ_grid).h5")
    mkpath(outdir)
    println(stdout)
    println(stderr)
    @info "Writing results" outfile
    flush(stdout)
    flush(stderr)
    HDF5.h5open(outfile, "w") do ff
        ff["rho"] = ρ_grid
        ff["Npoints"] = Np_actual
        ff["times_min"] = times_min
    end
    nothing
end

outdir = "results.CPU"
backend_short = CellListsBackend(CPU(), 2)
backend_long = NonuniformFFTsBackend(CPU(); σ = 1.5, m = HalfSupport(4))
Np_wanted_all = logrange(1e3, 1e6; length = 10)
for ρ_grid in (0.1,)
    run_benchmark_series(Np_wanted_all; ρ_grid, outdir, backend_short, backend_long)
end

# On CUDA (needs `using CUDA`):
if @isdefined(CUDABackend)
    outdir = "results.CUDABackend"
    backend_short = CellListsBackend(CUDABackend(), 2)
    backend_long = NonuniformFFTsBackend(CUDABackend(); σ = 1.5, m = HalfSupport(4))
    Np_wanted_all = logrange(1e4, 1e7; length = 10)
    for ρ_grid in (0.05, 0.1, 0.2)
        run_benchmark_series(Np_wanted_all; ρ_grid, outdir, backend_short, backend_long)
    end
end

# On CUDA with 2 GPUs:
if @isdefined(CUDABackend) && length(CUDA.devices()) >= 2
    outdir = "results.CUDABackend_2gpu"
    backend_short = CellListsBackend(CUDABackend(), 2; device = 2)
    backend_long = NonuniformFFTsBackend(CUDABackend(); device = 1, σ = 1.5, m = HalfSupport(4))
    Np_wanted_all = logrange(1e4, 1e7; length = 10)
    for ρ_grid in (0.05, 0.1, 0.2)
        run_benchmark_series(Np_wanted_all; ρ_grid, outdir, backend_short, backend_long)
    end
end

# On AMDGPU (needs `using AMDGPU`)
if @isdefined(ROCBackend)
    outdir = "results.ROCBackend"
    backend_short = CellListsBackend(ROCBackend(), 2)
    backend_long = NonuniformFFTsBackend(ROCBackend(); σ = 1.5, m = HalfSupport(4))
    Np_wanted_all = logrange(1e4, 1e7; length = 10)
    for ρ_grid in (0.05, 0.1, 0.2)
        run_benchmark_series(Np_wanted_all; ρ_grid, outdir, backend_short, backend_long)
    end
end

# On AMDGPU with 2 GPUs:
if @isdefined(ROCBackend) && length(AMDGPU.devices()) >= 2
    outdir = "results.ROCBackend_2gpu"
    backend_short = CellListsBackend(ROCBackend(), 2; device = 2)
    backend_long = NonuniformFFTsBackend(ROCBackend(); device = 1, σ = 1.5, m = HalfSupport(4))
    Np_wanted_all = logrange(1e4, 1e7; length = 10)
    for ρ_grid in (0.05, 0.1, 0.2)
        run_benchmark_series(Np_wanted_all; ρ_grid, outdir, backend_short, backend_long)
    end
end
