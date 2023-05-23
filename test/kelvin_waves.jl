# This tests the time evolution of a Kelvin wave.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm, normalize
using Optim: Optim
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart
using VortexPasta.Timestepping

# Initialise nearly straight vortex line with sinusoidal perturbation.
function init_vortex_line(; x, y, Lz = 2π, sign, A = 0.01, k::Int = 1,)
    tlims = (-0.5, 0.5)
    S(t) = SVector(
        x + A * sinpi(4π * k * t / Lz),
        y,
        (0.5 + sign * t) * Lz,
    )
    offset = setindex(zero(S(0.0)), sign * Lz, 3)
    (; x, y, Lz, sign, A, k, tlims, S, offset,)
end

function test_kelvin_waves(; Lz = 2π, A = 0.01, k = 1,)
    Lx = Ly = Lz
    lines = [
        init_vortex_line(; x = 0.25Lx, y = 0.25Ly, Lz, sign = +1, A = +A, k,),
        init_vortex_line(; x = 0.25Lx, y = 0.75Ly, Lz, sign = -1, A = -A, k,),
        init_vortex_line(; x = 0.75Lx, y = 0.25Ly, Lz, sign = -1, A = -A, k,),
        init_vortex_line(; x = 0.75Lx, y = 0.75Ly, Lz, sign = +1, A = +A, k,),
    ]

    N = 42
    filaments = map(lines) do line
        (; tlims, S, offset,) = line
        dζ = (tlims[2] - tlims[1]) / N
        ζs = range(tlims[1] + dζ / 2, tlims[2]; step = dζ)
        @assert length(ζs) == N && last(ζs) ≈ tlims[2] - dζ / 2
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod(); offset)
        # Filaments.init(ClosedFilament, S.(ζs), FiniteDiffMethod(2); offset)
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
        Ns = (1, 1, 1) .* 32
        kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
        α = kmax / 5
        rcut = 4 / α
        ParamsBiotSavart(;
            Γ, α, a, Δ, rcut, Ls, Ns,
            backend_short = NaiveShortRangeBackend(),
            backend_long = FINUFFTBackend(),
            quadrature_short = GaussLegendreQuadrature(4),
            quadrature_long = GaussLegendreQuadrature(4),
        )
    end

    # Determine timestep
    ℓ_min = minimum(filaments) do f
        local ts = knots(f)
        minimum(eachindex(segments(f))) do i
            ts[i + 1] - ts[i]
        end
    end
    β = Γ / 4π * (log(2 * ℓ_min / a) - Δ)
    dt_base = 1.8 * ℓ_min^2 / β        # this seems to give stability with RK4 (fails with 1.9)
    dt = T_kw / round(T_kw / dt_base)  # make sure that we exactly sample t = T_kw

    fs = copy.(filaments)
    tspan = (0.0, 2.2 * T_kw)

    jprobe = 5
    X_probe = eltype(fs[1])[]
    times = Float64[]

    function callback(iter)
        push!(times, iter.t)
        push!(X_probe, iter.fs[1][jprobe])
    end

    prob = @inferred VortexFilamentProblem(fs, tspan, params_bs)

    iter = @inferred init(
        prob, RK4();
        dt,
        adaptive = false,
        refinement = NoRefinement(),  # make sure that nodes don't "move" vertically due to refinement
        callback,
    )

    (; fs, vs,) = iter  # `vs` already contains the initial velocities

    # Analyse just one of the filaments.
    line = lines[1]
    fil = fs[1]
    vel = vs[1]

    @testset "Initial condition" begin
        @test iter.t == 0
        @test iter.nstep == 0
        sign_kw = -line.sign  # KWs rotate in the direction opposite to the vortex circulation
        for (X, v) ∈ zip(nodes(fil), vel)
            @test norm(v) ≈ abs(v[2])  # initial velocity is in the `y` direction
            δx = X[1] - line.x         # local perturbation
            v_expected = sign_kw * ω_kw * δx
            @test isapprox(v_expected, v[2]; rtol = 0.01)  # `rtol` can be further decreased if one increases the resolution `N` (things seem to converge)
        end
    end

    @info "Solving..."
    @time solve!(iter)

    # Analyse the trajectory of a single node
    @testset "Node trajectory" begin
        @test length(X_probe) == length(times) == iter.nstep + 1  # data include step 0
        # Analyse probed data
        xs = getindex.(X_probe, 1)
        ys = getindex.(X_probe, 2)
        zs = getindex.(X_probe, 3)

        # Position along the filament direction
        # @show abs.(zs .- zs[begin]) ./ zs[begin]
        @test all(z -> isapprox(z, zs[begin]; rtol = 1e-6), zs)  # the chosen node practically doesn't move vertically

        # Set origin at equilibrium position
        @. xs -= line.x
        @. ys -= line.y

        xmodel(t, p) = p[1] .* cos.(p[2] * t)
        ymodel(t, p) = p[1] .* sin.(p[2] * t)
        p0 = [line.A * 2, 3.2 * ω_kw]  # initial guess for parameters

        # Try to fit sinusoidal trajectory using Optim.
        # We minimise the square distance between data and model.
        # The particle swarm algorithm is much more robust for finding global
        # minima compared to gradient-based methods (such as LBFGS).
        xfit = Optim.optimize(p0, Optim.ParticleSwarm()) do p
            sum(zip(xs, times)) do (x, t)
                abs2(x - xmodel(t, p))
            end
        end
        yfit = Optim.optimize(p0, Optim.ParticleSwarm()) do p
            sum(zip(ys, times)) do (y, t)
                abs2(y - ymodel(t, p))
            end
        end

        # Verify that we found a good global minimum.
        # For some reason, the `y` fit is much better than the `x` one (which is still very good).
        Nt = length(times)
        @test xfit.minimum / abs(line.A * Nt) < 1e-9
        @test yfit.minimum / abs(line.A * Nt) < 1e-13

        # Compare amplitude of both signals
        Ax = abs(xfit.minimizer[1])
        Ay = abs(yfit.minimizer[1])
        ωx = abs(xfit.minimizer[2])
        ωy = abs(yfit.minimizer[2])

        @test isapprox(Ax, Ay; rtol = 1e-2)
        @test isapprox(ωx, ωy; rtol = 1e-4)

        # Check that we actually recover the KW frequency!
        @test isapprox(ωx, ω_kw; rtol = 2e-3)
        @test isapprox(ωy, ω_kw; rtol = 2e-3)
    end

    nothing
end

@testset "Kelvin waves" begin
    test_kelvin_waves()
end
