using Test
using Statistics: mean, std
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using UnicodePlots: lineplot
using LinearAlgebra: ⋅

@testset "Perturbed vortex ring" begin
    R = π / 4
    N = 32
    r_perturb(t) = 0.1 * sinpi(6t)  # mode m = 3
    S = @inferred define_curve(Ring(r = r_perturb); scale = R, translate = (π, π, 0))
    f = @inferred Filaments.init(S, ClosedFilament, N, CubicSplineMethod())

    Ls = (1, 1, 1) .* 2π
    Ns = (1, 1, 1) .* 32
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    α = kmax / 5
    rcut = 4 * sqrt(2) / α
    Γ = 1.2
    a = 1e-6
    Δ = 1/2
    params_bs = @inferred ParamsBiotSavart(;
        Γ, a, Δ,
        α, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = FINUFFTBackend(),
        quadrature_short = GaussLegendre(4),
        quadrature_long = GaussLegendre(4),
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

    times::Vector{Float64} = Float64[]
    energy_time::Vector{Float64} = Float64[]

    function callback(iter)
        (; nstep, t,) = iter.time
        if nstep == 0
            empty!(times)
            empty!(energy_time)
        end
        E = Diagnostics.kinetic_energy_from_streamfunction(iter; quad = GaussLegendre(4))
        # @show nstep, t, E
        push!(times, t)
        push!(energy_time, E)
        if @isdefined(Makie)
            fobs[] = iter.fs[1]  # update filament positions
            yield()
        end
        local (; fs, ψs, vs,) = iter
        if @isdefined(write_vtkhdf)  # requires loading VortexPasta.FilamentIO
            write_vtkhdf("ring_perturbed_$nstep.hdf", fs; refinement = 4) do io
                io["velocity"] = vs
                io["streamfunction"] = ψs
                ψt = similar(ψs, Float64)
                for (f, ψs, ψt) ∈ zip(fs, ψs, ψt)
                    for i ∈ eachindex(f, ψs, ψt)
                        ψt[i] = f[i, UnitTangent()] ⋅ ψs[i]
                    end
                    ψt[end + 1] = ψt[begin]  # assume it's a PaddedVector
                end
                io["streamfunction_t"] = ψt
            end
        end
        nothing
    end

    tspan = (0.0, 15 * R^2 / Γ)
    prob = @inferred VortexFilamentProblem([f], tspan, params_bs)
    iter = init(prob, SanduMRI33a(RK4(), 2); dt = 0.10, callback)
    @time solve!(iter)

    Emean = mean(energy_time)
    Estd = std(energy_time)
    @show Estd / Emean
    @test Estd / Emean < 4e-4

    if @isdefined(Makie)
        lines(times, energy_time)
    else
        plt = lineplot(
            times, energy_time;
            xlabel = "Time", ylabel = "Energy",
            title = "Perturbed vortex ring",
        )
        println(plt)
    end
end
