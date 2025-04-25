# Test restarts, based on trefoil knot reconnection.

using Test
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics

function test_restarts(::Type{T}; scheme = RK4(), dt_factor = T(0.5),) where {T}
    S = define_curve(TrefoilKnot(); translate = T(π), scale = T(π / 4))
    N = 64
    method = QuinticSplineMethod()
    f = Filaments.init(S, ClosedFilament{T}, N, method)
    fs_init = [f]

    params_bs = let L = T(2π)
        Ls = (1, 1, 1) .* L
        β = T(3.0)  # accuracy coefficient
        rcut = L / 3
        α = β / rcut
        kmax = 2α * β
        M = ceil(Int, kmax * L / π)
        Ns = (1, 1, 1) .* M
        ParamsBiotSavart(
            T;
            Γ = T(2.0), α, a = T(1e-6), Δ = T(1/4), rcut, Ls, Ns,
            backend_short = CellListsBackend(2),
            backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(3)),
            quadrature = GaussLegendre(3),
            lia_segment_fraction = T(0.2),
        )
    end

    times = T[]
    energy = T[]
    num_reconnections = Int[]
    num_rejections = Int[]

    function callback(iter)
        push!(times, iter.t)
        push!(energy, Diagnostics.kinetic_energy(iter))
        push!(num_reconnections, iter.stats.reconnection_count)
        push!(num_rejections, iter.nrejected)
        # @show iter.t, iter.stats.reconnection_count, iter.nrejected
        nothing
    end

    tsim = T(2.0)
    prob = @inferred VortexFilamentProblem(fs_init, tsim, params_bs)
    δ = Filaments.minimum_node_distance(prob.fs)
    d_crit = T(0.75) * δ
    refinement = RefineBasedOnSegmentLength(d_crit)
    reconnect = ReconnectBasedOnDistance(d_crit)
    dt = BiotSavart.kelvin_wave_period(params_bs, d_crit) * T(dt_factor)
    adaptivity = AdaptBasedOnVelocity(d_crit / 8; safety_factor = T(1.0))  # use large safety_factor to ensure rejections (usually bad!)

    iter = @inferred init(prob, scheme; dt, refinement, reconnect, adaptivity, callback)
    solve!(iter)

    N = lastindex(times)  # end of first simulation in output vectors

    # Now perform the same simulation with a checkpoint in-between.
    # We try to put the checkpoint after the first reconnections, when some timestep
    # rejections from the adaptivity criterion have already happened (these have some
    # memory that should be captured by the checkpoint).
    # Note that we will append values to existent times, energy, ... vectors.
    n = searchsortedlast(times, T(1.75))
    tsim_1 = times[n]
    @test num_reconnections[n] > 2
    @test num_rejections[n] > 10

    prob_1 = @inferred VortexFilamentProblem(fs_init, tsim_1, params_bs)
    iter_1 = @inferred init(prob_1, scheme; dt, refinement, reconnect, adaptivity, callback)
    solve!(iter_1)
    save_checkpoint("trefoil_checkpoint.vtkhdf", iter_1)

    @test iter_1.t == times[end]

    tsim_2 = tsim - iter_1.t
    checkpoint = @inferred load_checkpoint("trefoil_checkpoint.vtkhdf", T, method)
    prob_2 = @inferred VortexFilamentProblem(checkpoint, tsim_2, params_bs)
    @test prob_2.tspan[1] == iter_1.t
    iter_2 = @inferred init(prob_2, scheme; dt, refinement, reconnect, adaptivity, callback)
    solve!(iter_2)

    @test times[N] == times[end]
    @test energy[N] ≈ energy[end] rtol=1e-4  # note: the variability comes from threading; with 1 thread these are basically equal

    nothing
end

@testset "Checkpoints" begin
    test_restarts(Float64; scheme = RK4(), dt_factor = 0.5)
end
