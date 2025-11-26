# Define BenchmarkTools suite for CI
# https://github.com/MilesCranmer/AirspeedVelocity.jl

using VortexPasta
using VortexPasta.PredefinedCurves
using VortexPasta.Containers: VectorOfVectors
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Reconnections
using VortexPasta.Timestepping
using VortexPasta.Forcing
using Rotations: Rotations
using StableRNGs: StableRNG
using BenchmarkTools

using ThreadPinning
pinthreads(:cores)
threadinfo()

function generate_biot_savart_parameters(;
        Ls::NTuple{3, T}, β = 3.5,
        Ns::Union{Nothing, Dims{3}} = nothing,
        rcut::Union{Nothing, Real} = Ns === nothing ? min(Ls...) / 3 : nothing,
        kws...,
    ) where {T}
    Γ = 1.0
    a = 1e-7
    Δ = 1/2
    if Ns === nothing
        @assert rcut isa Real
        α = β / rcut
        kmax = 2α * β
        Ns = ceil.(Int, (kmax / π) .* Ls) .+ 2
    else
        @assert rcut === nothing
        ks_max = @. π * (Ns - 2) / Ls
        kmax = min(ks_max...)
        α = kmax / (2 * β)
        rcut = β / α
    end
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = NonuniformFFTsBackend(CPU(); σ = T(1.5), m = HalfSupport(4)),
        # backend_long = NonuniformFFTsBackend(CUDABackend(); σ = T(1.5), m = HalfSupport(4)),
        quadrature = GaussLegendre(3),
        quadrature_near_singularity = GaussLegendre(3),
        lia_segment_fraction = 0.2,
        kws...,
    )
end

abstract type VFMInitialCondition end

@kwdef struct SingleRing{T} <: VFMInitialCondition
    R::T
end

function generate_filaments(ring::SingleRing; Ls, l_res, method)
    (; R,) = ring
    T = eltype(Ls)
    N = ceil(Int, 2π * R / l_res)
    S = define_curve(Ring(); translate = Ls ./ 2, scale = R)
    f = Filaments.init(S, ClosedFilament{T}, N, method)
    [f]
end

@kwdef struct RandomRings{T} <: VFMInitialCondition
    nrings::Int
    R::T
end

function generate_filaments(rings::RandomRings; Ls, l_res, method)
    (; nrings, R,) = rings
    rng = StableRNG(42)
    N = ceil(Int, 2π * R / l_res)
    map(1:nrings) do _
        T = eltype(Ls)
        translate = Ls .* rand(rng, Vec3{T})
        rotate = rand(rng, Rotations.QuatRotation)
        S = define_curve(Ring(); scale = R, translate, rotate)
        Filaments.init(S, ClosedFilament{T}, N, method)
    end
end

## Create benchmark suite for AirspeedVelocity.jl (benchmarks on github CI)
const SUITE = BenchmarkGroup()

## Cell lists benchmarks
include(joinpath(@__DIR__, "cell_lists.jl"))
SUITE["CellLists"] = CellListsBenchmarks.main()
# results = run(SUITE["CellLists"])

## Prepare vortex configuration

Ls = (2π, 2π, 2π)
l_res = Ls[1] / 128  # typical line resolution
method = QuinticSplineMethod()
initial_condition = RandomRings(nrings = 500, R = π / 3)
fs = VectorOfVectors(generate_filaments(initial_condition; Ls, l_res, method))

## Define Biot–Savart benchmarks

Ns = (128, 128, 128)
params = generate_biot_savart_parameters(; Ls, Ns)
cache = BiotSavart.init_cache(params)
vs = similar(fs)
ψs = similar(fs)
fields = (velocity = vs, streamfunction = ψs)
BiotSavart.velocity_on_nodes!(vs, cache, fs)
BiotSavart.compute_on_nodes!(fields, cache, fs)
reset_timer!(cache.to)

SUITE["BiotSavart"]["velocity"] = @benchmarkable BiotSavart.velocity_on_nodes!($vs, $cache, $fs)
SUITE["BiotSavart"]["velocity + streamfunction"] = @benchmarkable BiotSavart.compute_on_nodes!($fields, $cache, $fs)

## Define reconnection benchmarks

let reconnect = ReconnectFast(l_res; max_passes = 10)
    cache_rec = Reconnections.init_cache(reconnect, fs, Ls)
    SUITE["Reconnections"]["ReconnectFast"] = @benchmarkable Reconnections.reconnect!($cache_rec, fs_rec) setup=(fs_rec = copy(fs))
end
let reconnect = ReconnectBasedOnDistance(l_res; max_passes = 10)
    cache_rec = Reconnections.init_cache(reconnect, fs, Ls)
    SUITE["Reconnections"]["ReconnectBasedOnDistance"] = @benchmarkable Reconnections.reconnect!($cache_rec, fs_rec) setup=(fs_rec = copy(fs))
end

## Define timestepping benchmarks

prob = VortexFilamentProblem(fs, 0.1, params)

l_min = 0.75 * l_res
adaptivity = AdaptBasedOnSegmentLength(0.5) | AdaptBasedOnVelocity(0.5 * l_min)
refinement = RefineBasedOnSegmentLength(l_min)
reconnect = ReconnectFast(l_min; max_passes = 10)
forcing = FourierBandForcingBS(; kmin = 0.1, kmax = 2.5, ε_target = 100.0, modify_length = false)
iter = init(prob, RK4(); dt = 0.01, adaptivity, refinement, reconnect, forcing)
step!(iter)
reset_timer!(iter.to)

SUITE["Timestepping"]["step!"] = @benchmarkable step!($iter) seconds=1000 evals=1 samples=10   # perform 10 timesteps
SUITE["Timestepping"]["forcing"] = @benchmarkable Timestepping.apply_forcing!(fields, $iter, fs, t, to) setup = begin
    fields = (; velocity = $iter.vL)
    fs = $iter.fs
    to = $iter.to
    t = $iter.t;
end

# results_timestepping = run(SUITE["Timestepping"])

## Example interactive usage:

# @time tune!(SUITE)  # not sure this is really needed (and takes a lot of time)
# results = run(SUITE)  # run full suite
# results_bs = run(SUITE["BiotSavart"])  # run BiotSavart benchmarks
# results_rec = run(SUITE["Reconnections"])
