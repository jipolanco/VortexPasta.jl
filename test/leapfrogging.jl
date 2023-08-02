# Test the time evolution of a system of two leapfrogging vortex rings.

using Test
using Statistics: mean, std
using Random
using LinearAlgebra: norm, normalize, ⋅
using UnicodePlots: lineplot
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Timestepping: VortexFilamentSolver

# Note: this is total energy within a unit cell.
# One should normalise by the cell volume to get energy per unit mass [L²T⁻²].
function kinetic_energy_from_streamfunction(
        ψs_all::AbstractVector, fs::AbstractVector{<:AbstractFilament},
        Γ::Real;
        quad = nothing,
    )
    prefactor = Γ / 2
    E = zero(prefactor)
    interpolate_tangent_component_only = true
    for (f, ψs) ∈ zip(fs, ψs_all)
        ts = knots(f)
        if quad === nothing
            for i ∈ eachindex(segments(f))
                ψ⃗ = ψs[i]
                s⃗′ = f[i, Derivative(1)]
                δt = (ts[i + 1] - ts[i - 1]) / 2
                E += (ψ⃗ ⋅ s⃗′) * δt
            end
        else
            Xoff = Filaments.end_to_end_offset(f)
            ψ_int = Filaments.change_offset(similar(f), zero(Xoff))
            if interpolate_tangent_component_only
                for i ∈ eachindex(ψs)
                    ψ_t = ψs[i] ⋅ f[i, UnitTangent()]
                    ψ_int[i] = (ψ_t, 0, 0)  # we only care about the first component
                end
            else
                copy!(nodes(ψ_int), ψs)
            end
            update_coefficients!(ψ_int; knots = knots(f))
            for i ∈ eachindex(segments(f))
                E += integrate(f, i, quad) do ζ
                    ψ⃗ = ψ_int(i, ζ)
                    if interpolate_tangent_component_only
                        ψ⃗[1]
                    else
                        s⃗′ = f(i, ζ, Derivative(1))
                        ψ⃗ ⋅ s⃗′
                    end
                end
            end
        end
    end
    E * prefactor
end

function kinetic_energy_from_streamfunction(iter::VortexFilamentSolver)
    (; fs, ψs,) = iter
    (; Γ,) = iter.prob.p.common
    quad = GaussLegendre(4)  # this doesn't seem to change much the results...
    kinetic_energy_from_streamfunction(ψs, fs, Γ; quad)
end

function total_vortex_length(fs)
    quad = GaussLegendre(4)
    L = 0.0
    for f ∈ fs
        for seg ∈ segments(f)
            L += integrate(seg, quad) do ζ
                norm(seg(ζ, Derivative(1)))
            end
        end
    end
    L
end

@testset "Leapfrogging vortex rings" begin
    ##
    # Grid-related parameters
    Ls = (1, 1, 1) .* 2π
    Ns = (1, 1, 1) .* 32
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    α = kmax / 5
    rcut = 4 * sqrt(2) / α

    # Physical vortex parameters
    Γ = 1.2
    a = 1e-6
    Δ = 1/4  # full core
    params_bs = @inferred ParamsBiotSavart(;
        Γ, a, Δ,
        α, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = FINUFFTBackend(),
        quadrature_short = GaussLegendre(4),
        quadrature_long = GaussLegendre(4),
    )

    # Initialise filaments
    N = 32
    R_init = π / 3
    zs_init = [0.9, 1.1] .* π
    rng = MersenneTwister(42)
    fs_init = map(zs_init) do z
        S = define_curve(Ring(); translate = (π / 20, π, z), scale = R_init)
        τs = collect(range(0, 1; length = 2N + 1)[2:2:2N])
        τs .+= rand(rng, N) ./ 3N  # slightly randomise locations
        Filaments.init(S, ClosedFilament, τs, CubicSplineMethod())
    end

    ##

    # Define function to be run at each simulation timestep
    times::Vector{Float64} = Float64[]
    energy_time::Vector{Float64} = Float64[]
    line_length::Vector{Float64} = Float64[]
    output_vtkhdf::Bool = false

    function callback(iter)
        (; fs, vs, ψs,) = iter
        (; t, dt, nstep,) = iter.time
        push!(times, t)

        # This can be useful for visualisation
        if output_vtkhdf
            write_vtkhdf("leapfrogging_$nstep.hdf", fs) do io
                write_point_data(io, "Velocity", vs)
                write_point_data(io, "Streamfunction", ψs)
            end
        end

        E = kinetic_energy_from_streamfunction(iter)
        L = total_vortex_length(fs)

        # @show nstep, t, dt, E
        push!(energy_time, E)
        push!(line_length, L)
    end

    # Initialise simulation
    fs = copy.(fs_init);
    tmax = 5 * R_init^2 / Γ  # enough time for a couple of "jumps"
    tspan = (0.0, tmax)
    prob = @inferred VortexFilamentProblem(fs, tspan, params_bs);
    iter = @inferred init(
        prob, RK4();
        dt = 0.025,  # will be changed by the adaptivity
        adaptivity = AdaptBasedOnSegmentLength(1.5),
        refinement = NoRefinement(),
        callback,
    );

    # Run simulation
    @time solve!(iter);

    # Check that the callback is called at the initial time
    @test first(times) == first(tspan)

    @testset "Energy conservation" begin
        energy_initial = first(energy_time)  # initial energy
        energy_normalised = energy_time ./ energy_initial
        energy_mean = mean(energy_normalised)
        energy_std = std(energy_normalised)
        let plt = lineplot(times, energy_normalised; xlabel = "Time", ylabel = "Energy")
            println(plt)
        end
        let plt = lineplot(times, line_length; xlabel = "Time", ylabel = "Length")
            println(plt)
        end
        @show energy_initial
        @show extrema(energy_normalised)
        @show energy_std
        @show energy_mean - 1
        @test energy_std < 1e-5
        @test isapprox(energy_mean, 1; rtol = 1e-5)
    end

    ##
end
