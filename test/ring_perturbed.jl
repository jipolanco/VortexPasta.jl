using Test
using Statistics: mean, std
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using UnicodePlots: lineplot
using LinearAlgebra: ⋅

function test_perturbed_vortex_ring()
    R = π / 4
    N = 48
    r_perturb(t) = 0.10 * sinpi(6t)  # mode m = 3
    z_perturb(t) = 0.02 * sinpi(8t)  # mode m = 4
    p = @inferred Ring(r = r_perturb, z = z_perturb)
    S = @inferred define_curve(p; scale = R, translate = (π, π, 0))
    f = @inferred Filaments.init(S, ClosedFilament, N, CubicSplineMethod())
    Filaments.redistribute_nodes!(f)

    Ls = (1, 1, 1) .* 2π
    Ns = (1, 1, 1) .* 32
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    β = 3.0
    α = kmax / 2β
    rcut = β / α
    Γ = 1.2
    a = 1e-6
    Δ = 1/2
    params_bs = @inferred ParamsBiotSavart(;
        Γ, a, Δ,
        α, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = FINUFFTBackend(tol = 1e-6),
        quadrature = GaussLegendre(2),
    )

    if @isdefined(Makie)
        fig = Figure()
        ax = Axis3(fig[1, 1]; aspect = :data)
        hidespines!(ax)
        wireframe!(ax, Rect(0, 0, 0, Ls...); color = :grey, linewidth = 0.5)
        fobs = Observable(f)
        plot!(ax, fobs; refinement = 4, markersize = 4)
        display(fig)
    end

    times = Float64[]
    energy_time = Float64[]
    helicity_time = Float64[]
    stretching_time = Float64[]
    tsf = TimeSeriesFile()

    function callback(iter)
        (; nstep, t,) = iter.time
        if nstep == 0
            empty!(times)
            empty!(energy_time)
            empty!(helicity_time)
            empty!(tsf)
        end
        quad = GaussLegendre(3)
        E = Diagnostics.kinetic_energy_from_streamfunction(iter; quad)
        H = Diagnostics.helicity(iter; quad)
        dLdt = Diagnostics.stretching_rate(iter; quad)
        # @show nstep, t, E
        push!(times, t)
        push!(energy_time, E)
        push!(helicity_time, H)
        push!(stretching_time, dLdt)
        if @isdefined(Makie)
            fobs[] = iter.fs[1]  # update filament positions
            yield()
        end
        local (; fs, ψs, vs,) = iter
        if nstep % 100 == 0  # don't output very often; just for tests
            let fname = "ring_perturbed_$nstep.vtkhdf"
                write_vtkhdf(fname, fs; refinement = 4, periods = Ls) do io
                    io["velocity"] = vs
                    io["streamfunction"] = ψs
                    ψt = map(f -> similar(nodes(f), Float64), fs)
                    for (f, ψs, ψt) ∈ zip(fs, ψs, ψt)
                        for i ∈ eachindex(f, ψs, ψt)
                            ψt[i] = f[i, UnitTangent()] ⋅ ψs[i]
                        end
                        ψt[end + 1] = ψt[begin]  # assume it's a PaddedVector
                    end
                    io["streamfunction_t"] = ψt
                end
                tsf[t] = fname
            end
        end
        nothing
    end

    tspan = (0.0, 15 * R^2 / Γ)
    prob = @inferred VortexFilamentProblem([f], tspan, params_bs)

    # When we only use the AdaptBasedOnVelocity criterion, then `init` will complete it with
    # a `MaximumTimestep` criterion. We check this below.
    adaptivity_in = AdaptBasedOnVelocity(100.0)  # huge value of δ to make sure it doesn't actually limit the timestep
    d_min = Filaments.minimum_node_distance(prob.fs)
    dt_kw = BiotSavart.kelvin_wave_period(params_bs, d_min)
    dt = 10 * dt_kw
    iter = @inferred init(prob, SanduMRI33a(RK4(), 3); dt, adaptivity = adaptivity_in, callback)

    let adaptivity = iter.adaptivity
        @test adaptivity isa Timestepping.CombinedAdaptivityCriteria
        (; criteria,) = adaptivity
        criteria :: Tuple{Vararg{Timestepping.AdaptivityCriterion}}  # `criteria` is a tuple of criteria
        @test length(criteria) == 2
        @test criteria[1] === adaptivity_in
        @test criteria[2] === MaximumTimestep(dt)
    end

    @time solve!(iter)
    save("ring_perturbed.vtkhdf.series", tsf)



    @test iter.dt == dt  # dt was kept to the initial value (because AdaptBasedOnVelocity criterion gives larger dt)

    # Check that velocity and streamfunction can be interpolated, and that the interpolated
    # values on the nodes match the actual values at the nodes.
    @testset "Interpolate velocity & streamfunction" begin
        let fs = iter.fs, i = lastindex(fs), f = fs[i]
            for us ∈ (iter.vs[i], iter.ψs[i])
                @test length(us) == length(f)
                @test Filaments.knots(us) == Filaments.knots(f)
                @test Filaments.parametrisation(us) === Filaments.parametrisation(f)
                @test us[4] ≈ us(4, 0.0)
                @test us[5] ≈ us(4, 1.0)
            end
        end
    end

    # Helicity should be zero at all times.
    Hnorm = helicity_time ./ Γ^2
    # @show extrema(Hnorm)
    @test maximum(abs, Hnorm) < 2e-5

    Emean = mean(energy_time)
    Estd = std(energy_time)
    @show Estd / Emean
    @test Estd / Emean < 1e-4

    if @isdefined(Makie)
        lines(times, energy_time)
        display(current_figure())
    else
        plt = lineplot(
            times, energy_time;
            xlabel = "Time", ylabel = "Energy",
            title = "Perturbed vortex ring",
        )
        println(plt)
    end

    nothing
end

@testset "Perturbed vortex ring" begin
    test_perturbed_vortex_ring()
end
