# This tests the time evolution of a Kelvin wave.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm, normalize, ⋅
using UnicodePlots: lineplot
using Optim: Optim
using VortexPasta.PredefinedCurves: define_curve, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics

using JET: @test_opt
using FINUFFT: FINUFFT  # for JET only

# Initialise nearly straight vortex line with sinusoidal perturbation.
function init_vortex_line(; x, y, Lz = 2π, sign, A = 0.01, k::Int = 1,)
    tlims = (0, 1)
    xfun(t) = A * sinpi(2 * k * t)
    p = PeriodicLine(; x = xfun,)
    S = define_curve(
        p;
        scale = (1, 1, Lz),  # make the line 2π-periodic in Z
        translate = (x, y, Lz / 2),
        orientation = sign,
    )
    offset = S(1) - S(0)
    @assert offset ≈ SVector(0, 0, sign * Lz)
    (; x, y, Lz, sign, A, k, tlims, S, offset,)
end

# Explicit RK methods
dt_factor(::RK4) = 1.5    # this factor seems to give stability with RK4
dt_factor(::DP5) = 1.3    # I'd expect DP5 to allow a larger timestep than RK4, but that doesn't seem to be the case...
dt_factor(::SSPRK33) = 1.1
dt_factor(::Euler) = 0.08  # Euler needs a really small dt to stay stable, and accuracy is quite bad!!
dt_factor(::Midpoint) = 0.4

# IMEX-RK methods
dt_factor(::IMEXEuler) = 0.25
dt_factor(::Ascher343) = 1.2
dt_factor(::KenCarp3) = 1.2
dt_factor(::KenCarp4) = 1.5  # we could increase this but then it takes more iterations and time to solve the implicit problem!

dt_factor(::MultirateMidpoint) = 4.0
dt_factor(::SanduMRI33a) = 5.0
dt_factor(::SanduMRI45a) = 8.0

function test_kelvin_waves(
        scheme = RK4();
        method = QuinticSplineMethod(), Lz = 2π, A = 0.01, k = 1,
        quad = GaussLegendre(4),
        regularise_binormal = Val(false),  # enables larger timesteps (and accuracy loss...)
    )
    test_jet = true
    Lx = Ly = Lz
    lines = [
        init_vortex_line(; x = 0.25Lx, y = 0.25Ly, Lz, sign = +1, A = +A, k,),
        init_vortex_line(; x = 0.25Lx, y = 0.75Ly, Lz, sign = -1, A = -A, k,),
        init_vortex_line(; x = 0.75Lx, y = 0.25Ly, Lz, sign = -1, A = -A, k,),
        init_vortex_line(; x = 0.75Lx, y = 0.75Ly, Lz, sign = +1, A = +A, k,),
    ]

    N = 16
    filaments = map(lines) do line
        Filaments.init(line.S, ClosedFilament, N, method)
    end

    Γ = 2.0
    a = 1e-6
    Δ = 1/4  # full core (same as Schwarz 1985)

    # Expected KW frequency, from Schwarz 1985 (eq. A14)
    # In fact, this formulation (with 1/2 - Δ) is supposed to fit both Δ = 1/4
    # (full core, e.g. Schwarz 1985) and Δ = 1/2 (hollow core, e.g. Donnelly 1991) cases.
    γ = MathConstants.eulergamma  # Euler–Mascheroni constant
    ω_kw = Γ * k^2 / (4 * π) * (
        log(2 / (k * a)) - γ + 1/2 - Δ
    )
    T_kw = 2π / ω_kw

    params_bs = let
        @assert Lx == Ly == Lz
        Ls = (Lx, Ly, Lz)
        Ns = (1, 1, 1) .* 16
        kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
        α = kmax / 5
        rcut = 5 / α
        ParamsBiotSavart(;
            Γ, α, a, Δ, rcut, Ls, Ns,
            # backend_short = CellListsBackend(2),
            backend_short = NaiveShortRangeBackend(),
            backend_long = FINUFFTBackend(tol = 1e-6),
            quadrature = quad,
            regularise_binormal,
            # For NoQuadrature, it doesn't make much sense to use lia_segment_fraction ≠ 1.
            lia_segment_fraction = quad === NoQuadrature() ? 1.0 : 0.1,
        )
    end

    fs = copy.(filaments)
    tmax = if scheme isa Euler
        1.0  # explicit Euler needs a very small dt, hence we reduce the time...
    else
        1.5
    end
    tspan = (0.0, tmax * T_kw)

    jprobe = 5
    X_probe = eltype(fs[1])[]
    times = Float64[]
    energy_time = Float64[]

    function callback(iter)
        (; t, nstep,) = iter
        push!(times, t)
        push!(X_probe, iter.fs[1][jprobe])
        T = iter.prob.tspan[2]
        local (; fs, ψs,) = iter
        local (; Γ, Ls,) = iter.prob.p.common
        local quad = nothing  # this doesn't seem to change much the results...
        E = Diagnostics.kinetic_energy_from_streamfunction(ψs, fs, Γ, Ls; quad)
        # if nstep % 10 == 0
        #     @show nstep, t/T, E
        # end
        # write_vtkhdf("kw_$nstep.vtkhdf", fs; refinement = 4) do io
        #     io["Streamfunction"] = ψs
        # end
        push!(energy_time, E)
    end

    prob = @inferred VortexFilamentProblem(fs, tspan, params_bs)

    # This tests `show(::IO, prob)`.
    @test startswith(repr(prob), "VortexFilamentProblem with fields:\n")

    factor = dt_factor(scheme)
    if method isa CubicSplineMethod
        factor *= 0.8  # CubicSplineMethod requires slightly smaller timestep
    end

    iter = @inferred init(
        prob, scheme;
        dtmin = factor * T_kw * 1e-4,
        dt = 1.0,  # will be changed by the adaptivity
        fast_term = ShortRangeTerm(),  # allows to increase the timestep with splitting methods
        adaptivity = AdaptBasedOnSegmentLength(factor),
        refinement = NoRefinement(),  # make sure that nodes don't "move" vertically due to refinement
        callback,
    )

    # This tests `show(::IO, iter)`.
    @test startswith(repr(iter), "VortexFilamentSolver with fields:\n")
    @test startswith(repr(iter.time), "TimeInfo:\n")

    let
        ks, Ek = @inferred Diagnostics.energy_spectrum(iter)
        Δk = step(ks)
        E = Diagnostics.kinetic_energy_from_streamfunction(iter; quad = GaussLegendre(2))
        @test E/10 < sum(Ek) * Δk < E
        # plt = lineplot(ks[2:end], Ek[2:end]; xscale = log10, yscale = log10)
        # println(plt)
    end

    if test_jet
        @test_opt ignored_modules=(Base, FINUFFT) step!(iter)
    end

    (; fs, vs, ψs,) = iter  # `vs` already contains the initial velocities

    # Analyse just one of the filaments.
    line = lines[1]
    fil = fs[1]
    vel = vs[1]

    # This shouldn't depend on the scheme, but it can depend on the discretisation
    # parameters.
    @testset "Initial condition" begin
        @test iter.time.t == 0
        @test iter.time.nstep == 0
        sign_kw = -line.sign  # KWs rotate in the direction opposite to the vortex circulation
        for (X, v) ∈ zip(nodes(fil), vel)
            @test norm(v) ≈ abs(v[2])  # initial velocity is in the `y` direction
            δx = X[1] - line.x         # local perturbation
            v_expected = sign_kw * ω_kw * δx
            # `rtol` can be further decreased if one increases the resolution `N` (things seem to converge)
            # @show (v_expected - v[2]) / v_expected
            rtol = if method isa CubicSplineMethod
                0.03
            else
                0.01
            end
            @test isapprox(v_expected, v[2]; rtol)
        end
    end

    tlast = tspan[2]
    @info "Solving with $scheme..." dt_initial = iter.time.dt tlast/T_kw
    @time solve!(iter)
    println(iter.to)

    # Check that the callback is called at the initial time
    @test first(times) == first(tspan)

    let Enorm = energy_time ./ first(energy_time)
        plt = lineplot(
            times, Enorm;
            xlabel = "Time", ylabel = "Energy",
            title = "Kelvin waves / $scheme",
        )
        println(plt)
    end

    # Check energy conservation
    iseuler = scheme isa Euler || scheme isa IMEXEuler  # reduced precision of results
    @testset "Energy conservation" begin
        energy_initial = first(energy_time)  # initial energy
        energy_last = last(energy_time)
        energy_std = std(energy_time)
        @show energy_initial
        @show energy_std / energy_initial
        @show energy_last / energy_initial - 1
        if iseuler
            @test energy_std / energy_initial < 1e-6
            @test isapprox(energy_last, energy_initial; rtol = 1e-5)
        elseif method isa CubicSplineMethod || quad === NoQuadrature()
            @test energy_std / energy_initial < 2e-9
            @test isapprox(energy_last, energy_initial; rtol = 2e-9)
        else
            @test energy_std / energy_initial < 1e-9
            @test isapprox(energy_last, energy_initial; rtol = 2e-9)
        end
    end

    # Analyse the trajectory of a single node
    @testset "Node trajectory" begin
        @test length(X_probe) == length(times) == iter.time.nstep + 1  # data include step 0
        # Analyse probed data
        xs = getindex.(X_probe, 1)
        ys = getindex.(X_probe, 2)
        zs = getindex.(X_probe, 3)

        # Position along the filament direction
        # @show maximum(z -> abs(z - zs[begin]) / zs[begin], zs)
        let rtol = iseuler ? 1e-5 : 1e-6
            @test all(z -> isapprox(z, zs[begin]; rtol), zs)  # the chosen node practically doesn't move vertically
        end

        # Set origin at equilibrium position
        @. xs -= line.x
        @. ys -= line.y

        A = abs(line.A)
        @test std(xs) < A  # if this fails is likely because the solver diverged
        @test std(ys) < A

        xmodel(t, p) = p[1] * cos(p[2] * t)
        ymodel(t, p) = p[1] * sin(p[2] * t)

        # Initial guess for model parameters.
        # For gradient-based methods (such as LBFGS), the guessed frequency
        # must be quite close to the expected one, or otherwise the method
        # converges to a different local minimum.
        p0 = [2 * A, 1.1 * ω_kw]

        # Try to fit sinusoidal trajectory using Optim.
        # We minimise the square distance between data and model.
        xfit = Optim.optimize(p0, Optim.LBFGS(); autodiff = :forward) do p
            sum(zip(xs, times)) do (x, t)
                abs2(x - xmodel(t, p))
            end
        end
        yfit = Optim.optimize(p0, Optim.LBFGS(); autodiff = :forward) do p
            sum(zip(ys, times)) do (y, t)
                abs2(y - ymodel(t, p))
            end
        end

        # Verify that we found a good global minimum.
        A_tol = iseuler ? 2e-6 : 2e-13
        Nt = length(times)
        # @show xfit.minimum / (A * Nt)
        # @show yfit.minimum / (A * Nt)
        @test xfit.minimum / (A * Nt) < A_tol
        @test yfit.minimum / (A * Nt) < A_tol

        # Compare amplitude of both signals
        Ax = abs(xfit.minimizer[1])
        Ay = abs(yfit.minimizer[1])
        ωx = abs(xfit.minimizer[2])
        ωy = abs(yfit.minimizer[2])

        # @show abs(Ax - Ay) / Ay
        # @show abs(ωx - ωy) / ωy
        @test isapprox(Ax, Ay; rtol = 1e-2)
        @test isapprox(ωx, ωy; rtol = iseuler ? 2e-3 : 1e-6)

        # Check that we actually recover the KW frequency!
        # @show abs(ωy - ω_kw) / ω_kw
        rtol = method isa CubicSplineMethod ? 2e-2 : quad === NoQuadrature() ? 6e-3 : 1e-3
        @test isapprox(ωx, ω_kw; rtol)
        @test isapprox(ωy, ω_kw; rtol)
    end

    nothing
end

@testset "Kelvin waves" begin
    schemes = (
        RK4(), KenCarp3(),
        SanduMRI33a(RK4(), 1),
        # SanduMRI45a(RK4(), 1),
    )
    @testset "Scheme: $scheme" for scheme ∈ schemes
        test_kelvin_waves(scheme; method = QuinticSplineMethod())
    end
    @testset "NoQuadrature()" begin
        test_kelvin_waves(RK4(); quad = NoQuadrature())
    end
    @testset "CubicSplineMethod()" begin
        test_kelvin_waves(SanduMRI33a(RK4(), 1); method = CubicSplineMethod())
    end
    @testset "FourierMethod()" begin
        test_kelvin_waves(SanduMRI33a(RK4(), 1); method = FourierMethod())
    end
    @testset "FiniteDiffMethod()" begin
        test_kelvin_waves(RK4(); method = FiniteDiffMethod())
    end
end
