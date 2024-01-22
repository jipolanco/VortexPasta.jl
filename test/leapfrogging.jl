# Test the time evolution of a system of two leapfrogging vortex rings.

using Test
using Statistics: mean, std
using Random
using LinearAlgebra: norm, normalize, ⋅, ×
using UnicodePlots: lineplot
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Timestepping: VortexFilamentSolver
using VortexPasta.Diagnostics

using JET: JET
using FINUFFT: FINUFFT  # for JET only

function init_ring_filaments(R_init; method = CubicSplineMethod(), noise = 1/3)
    N = 32
    zs_init = [0.9, 1.1] .* π
    rng = MersenneTwister(42)
    map(zs_init) do z
        S = define_curve(Ring(); translate = (π / 20, π, z), scale = R_init)
        τs = collect(range(0, 1; length = 2N + 1)[2:2:2N])
        τs .+= rand(rng, N) .* (noise / N)  # slightly randomise locations
        Filaments.init(S, ClosedFilament, τs, method)
    end
end

function vortex_ring_squared_radius(f::AbstractFilament)
    quad = GaussLegendre(4)
    # Note: this is the impulse normalised by the vortex circulation Γ and the density ρ.
    # For a vortex ring on the XY plane, this should be equal to (0, 0, A) where A = πR² is
    # the vortex "area".
    i⃗ = Diagnostics.vortex_impulse(f; quad)
    A = i⃗[3]
    # A = norm(i⃗)
    A / π
end

# Timestep factors required for conservation tests to pass
dt_factor(::RK4) = 1.5
dt_factor(::DP5) = 1.5
dt_factor(::SSPRK33) = 1.0
dt_factor(::Euler) = 0.08
dt_factor(::Midpoint) = 0.4

dt_factor(::KenCarp3) = 1.3
dt_factor(::KenCarp4) = 2.3
dt_factor(::Ascher343) = 1.3

dt_factor(::MultirateMidpoint) = 0.8
dt_factor(::SanduMRI33a) = 4.0
dt_factor(::SanduMRI45a) = 4.0

# NOTE: this factor ensures stability, but *not* accuracy of IMEXEuler. The factor should
# actually be ~0.1 to get similar results to other schemes. With 0.6, one can clearly see
# the difference in results from the plots of the total vortex length.
dt_factor(::IMEXEuler) = 0.6

function test_leapfrogging_rings(
        prob, scheme;
        R_init,  # initial ring radii
        refinement,
        label = string(scheme),
    )
    test_jet = true

    # Define callback function to be run at each simulation timestep
    times = Float64[]
    energy_time = Float64[]
    impulse_time = Vec3{Float64}[]  # not used, just for inference tests (JET)
    line_length = Float64[]
    sum_of_squared_radii = Float64[]

    function callback(iter)
        local (; fs, vs, ψs, t, dt, nstep,) = iter
        push!(times, t)

        # This can be useful for visualisation
        # write_vtkhdf("leapfrogging_$nstep.vtkhdf", fs) do io
        #     io["Velocity"] = vs
        #     io["Streamfunction"] = ψs
        # end

        quad = GaussLegendre(4)
        E = Diagnostics.kinetic_energy_from_streamfunction(iter; quad)
        L = Diagnostics.filament_length(iter; quad)
        p⃗ = Diagnostics.vortex_impulse(iter; quad)

        # R²_all = @inferred sum(vortex_ring_squared_radius, fs)  # inference randomly fails on Julia 1.10-beta1...
        R²_all = 0.0
        for f ∈ fs
            R²_all += vortex_ring_squared_radius(f)
        end

        # @show nstep, t, dt, E
        push!(energy_time, E)
        push!(line_length, L)
        push!(impulse_time, p⃗)
        push!(sum_of_squared_radii, R²_all)
    end

    if test_jet
        JET.@test_opt ignored_modules=(Base, FINUFFT) init(prob, scheme; dt = 0.01)
        JET.@test_call ignored_modules=(Base, FINUFFT) init(prob, scheme; dt = 0.01)
    end

    l_min = minimum_knot_increment(prob.fs)

    # The 2π is simply because the factors were originally tuned with an old definition of
    # AdaptBasedOnSegmentLength.
    factor = dt_factor(scheme) * 2π
    if Filaments.discretisation_method(eltype(prob.fs)) isa QuinticSplineMethod
        # Quintic splines seem to need a smaller timestep...
        factor *= 0.3
    end
    adaptivity =
        AdaptBasedOnSegmentLength(factor) |
        AdaptBasedOnVelocity(10 * l_min) |  # usually inactive; just for testing
        AdaptBasedOnVelocity(20 * l_min)    # usually inactive; just for testing

    iter = @inferred init(
        prob, scheme;
        dt = 0.025,  # will be changed by the adaptivity
        dtmin = 0.005 * dt_factor(scheme),
        alias_u0 = false,  # don't overwrite fs_init
        adaptivity,
        refinement,
        callback,
    )

    if test_jet
        JET.@test_opt ignored_modules=(Base, FINUFFT) step!(iter)
        JET.@test_call ignored_modules=(Base, FINUFFT) step!(iter)
    end

    @info "Solving with $scheme..." dt_initial = iter.time.dt prob.tspan

    # Run simulation
    @time solve!(iter)

    println(iter.to)

    # Check that the callback is called at the initial time
    @test first(times) == first(prob.tspan)

    iseuler = scheme isa Euler || scheme isa IMEXEuler  # reduced precision of results
    @testset "Energy & impulse conservation" begin
        # With refinement we lose a tiny bit of precision (but still very acceptable!).
        rtol_energy = refinement === NoRefinement() ? 5e-5 : 1e-4
        rtol_impulse = 2e-5
        if iseuler
            rtol_energy *= 100
            rtol_impulse *= 200
        elseif scheme isa Midpoint
            rtol_energy *= 2
            rtol_impulse *= 2
        end

        energy_initial = first(energy_time)  # initial energy
        energy_normalised = energy_time ./ energy_initial

        # Impulse of a vortex ring is I = πR² × Γ
        # First, check that the squared radii are correctly estimated via the impulse.
        fs_init = prob.fs  # initial condition
        R²_sum_initial = length(fs_init) * R_init^2
        @test isapprox(R²_sum_initial, first(sum_of_squared_radii); rtol = 1e-5)

        impulse_normalised = sum_of_squared_radii ./ first(sum_of_squared_radii)

        let
            plt = lineplot(
                times, energy_normalised;
                xlabel = "Time", ylabel = "Energy",
                title = "Leapfrogging rings / " * label,
            )
            println(plt)
        end
        let
            plt = lineplot(times, line_length; xlabel = "Time", ylabel = "Length")
            println(plt)
        end
        let
            plt = lineplot(times, impulse_normalised; xlabel = "Time", ylabel = "Impulse")
            println(plt)
        end

        energy_mean = mean(energy_normalised)
        energy_std = std(energy_normalised)
        @show energy_initial
        @show extrema(energy_normalised)
        @show energy_std (energy_mean - 1)

        @test energy_std < rtol_energy
        @test isapprox(energy_mean, 1; rtol = rtol_energy)

        impulse_mean = mean(impulse_normalised)
        impulse_std = std(impulse_normalised)

        @show impulse_std (impulse_mean - 1)
        @test impulse_std < rtol_impulse
        @test isapprox(impulse_mean, 1; rtol = rtol_impulse)
    end

    nothing
end

##

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
        quadrature = GaussLegendre(4),
        regularise_binormal = Val(true),
    )

    # Check overloaded getproperty and propertynames for ParamsBiotSavart.
    @test α == @inferred (p -> p.α)(params_bs)
    @test :α ∈ @inferred propertynames(params_bs)

    # Initialise simulation
    R_init = π / 3
    fs_init = init_ring_filaments(R_init)
    tmax = 5 * R_init^2 / Γ  # enough time for a couple of "jumps"
    tspan = (0.0, tmax)
    prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs)

    @testset "VortexFilamentSolver" begin
        iter = @inferred init(prob, Midpoint(); dt = 0.01)
        @test iter isa VortexFilamentSolver

        # Check overloaded getproperty and propertynames for VortexFilamentSolver.
        @test iter.time.dt === @inferred (p -> p.dt)(iter)
        @test :dt ∈ @inferred propertynames(iter)
    end

    # Test diagnostics with/out quadratures here
    @testset "Diagnostics" begin
        rtol = 2e-3
        iter = init(prob, Midpoint(); dt = 0.01)
        @test isapprox(
            @inferred(Diagnostics.kinetic_energy_from_streamfunction(iter; quad = nothing)),
            @inferred(Diagnostics.kinetic_energy_from_streamfunction(iter; quad = GaussLegendre(2)));
            rtol,
        )
        @test isapprox(
            @inferred(Diagnostics.filament_length(iter; quad = nothing)),
            @inferred(Diagnostics.filament_length(iter; quad = GaussLegendre(2)));
            rtol,
        )
        @test isapprox(
            @inferred(Diagnostics.vortex_impulse(iter; quad = nothing)),
            @inferred(Diagnostics.vortex_impulse(iter; quad = GaussLegendre(2)));
            rtol,
        )
        let
            ks, Ek = @inferred Diagnostics.energy_spectrum(iter)
            Δk = step(ks)
            E = Diagnostics.kinetic_energy_from_streamfunction(iter; quad = GaussLegendre(2))
            # The energy spectrum doesn't include small-scale energy, so its integral
            # underestimates the total energy.
            @test E/10 < sum(Ek) * Δk < E
            # plt = lineplot(ks[2:end], Ek[2:end]; xscale = log10, yscale = log10)
            # println(plt)
            Ek_alt = similar(Ek)
            Diagnostics.energy_spectrum!(Ek_alt, ks, iter)
            @test Ek_alt == Ek
        end
    end

    δ = @inferred minimum_knot_increment(fs_init)
    @test δ ≈ @inferred minimum_node_distance(fs_init)  # almost always true, except for FourierMethod
    refinement = RefineBasedOnSegmentLength(0.99 * δ)
    schemes = (
        RK4(), DP5(), SSPRK33(),
        IMEXEuler(), KenCarp4(),
        KenCarp3(), Ascher343(),
        # Euler(),  # too slow!
        Midpoint(),
        MultirateMidpoint(32),
        SanduMRI33a(12),
        SanduMRI33a(CrankNicolson(), 2),
        SanduMRI45a(6),
    )

    ##

    @testset "Scheme: $scheme" for scheme ∈ schemes
        test_leapfrogging_rings(prob, scheme; R_init, refinement)
    end

    @testset "QuinticSplineMethod" begin
        local scheme = SanduMRI33a(RK4(), 2)
        local fs_init = @inferred init_ring_filaments(R_init; method = QuinticSplineMethod())
        local prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs)
        test_leapfrogging_rings(prob, scheme; R_init, refinement, label = "QuinticSplineMethod")
    end

    @testset "RK4 + no refinement" begin
        test_leapfrogging_rings(prob, RK4(); R_init, refinement = NoRefinement(), label = "NoRefinement")
    end

    @testset "FourierMethod (with noise = 0)" begin
        local scheme = SanduMRI33a(RK4(), 2)
        local fs_init = @inferred init_ring_filaments(R_init; method = FourierMethod(), noise = 0.0)
        local prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs)
        test_leapfrogging_rings(prob, scheme; R_init, refinement = NoRefinement(), label = "FourierMethod")
    end
end
