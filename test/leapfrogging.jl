# Test the time evolution of a system of two leapfrogging vortex rings.

using Test
using Statistics: mean, std
using LinearAlgebra: norm, normalize, ⋅
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Timestepping: VortexFilamentSolver

function kinetic_energy_from_streamfunction(
        ψs_all::AbstractVector, fs::AbstractVector{<:AbstractFilament},
        Γ::Real, Ls::NTuple{3};
        quad = nothing,
    )
    prefactor = Γ / (2 * prod(Ls))
    E = zero(prefactor)
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
            copy!(nodes(ψ_int), ψs)
            update_coefficients!(ψ_int; knots = knots(f))
            for i ∈ eachindex(segments(f))
                E += integrate(f, i, quad) do ζ
                    ψ⃗ = ψ_int(i, ζ)
                    s⃗′ = f(i, ζ, Derivative(1))
                    ψ⃗ ⋅ s⃗′
                end
            end
        end
    end
    E * prefactor
end

function kinetic_energy_from_streamfunction(iter::VortexFilamentSolver)
    (; fs, ψs,) = iter
    (; Γ, Ls,) = iter.prob.p.common
    quad = GaussLegendre(4)  # this doesn't seem to change much the results...
    kinetic_energy_from_streamfunction(ψs, fs, Γ, Ls; quad)
end

@testset "Leapfrogging vortex rings" begin
    # Grid-related parameters
    Ls = (1, 1, 1) .* 2π
    Ns = (1, 1, 1) .* 48
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    α = kmax / 6
    rcut = 4 * sqrt(2) / α
    @assert rcut < minimum(Ls) / 2  # required by the implementation

    # Physical vortex parameters
    Γ = 1.2
    a = 1e-4
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
    fs_init = map(zs_init) do z
        S = define_curve(Ring(); translate = (π, π, z), scale = R_init)
        Filaments.init(S, ClosedFilament, N, CubicSplineMethod())
    end

    ##

    # Define function to be run at each simulation timestep
    times::Vector{Float64} = Float64[]
    energy_time::Vector{Float64} = Float64[]
    output_vtkhdf::Bool = false
    function callback(iter)
        (; fs, vs, ψs,) = iter
        (; t, nstep,) = iter.time
        push!(times, t)

        # This can be useful for visualisation
        if output_vtkhdf
            write_vtkhdf("leapfrogging_$nstep.hdf", fs) do io
                write_point_data(io, "Velocity", vs)
                write_point_data(io, "Streamfunction", ψs)
            end
        end

        E = kinetic_energy_from_streamfunction(iter)
        push!(energy_time, E)
    end

    # Initialise simulation
    fs = copy.(fs_init);
    tmax = 4 * R_init^2 / Γ  # enough time for a couple of "jumps"
    tspan = (0.0, tmax)
    prob = @inferred VortexFilamentProblem(fs, tspan, params_bs);
    iter = @inferred init(
        prob, RK4();
        dt = 1.0,  # will be changed by the adaptivity
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
        energy_mean = mean(energy_time)
        energy_last = last(energy_time)
        energy_std = std(energy_time)
        @show energy_initial
        @show energy_std / energy_initial
        @show energy_last / energy_initial - 1
        @show energy_mean / energy_initial - 1
        @test energy_std / energy_initial < 1e-5
        @test isapprox(energy_mean, energy_initial; rtol = 1e-5)
        @test isapprox(energy_last, energy_initial; rtol = 1e-5)
    end
end
