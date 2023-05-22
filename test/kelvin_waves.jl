# This tests the time evolution of a Kelvin wave.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm
using FFTW: FFTW
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
    (; x, y, Lz, sign, tlims, S, offset,)
end

function test_kelvin_waves(; Lz = 2π, A = 0.01, k = 2,)
    Lx = Ly = Lz
    lines = [
        init_vortex_line(; x = 0.25Lx, y = 0.25Ly, Lz, sign = +1, A = +A, k,),
        init_vortex_line(; x = 0.25Lx, y = 0.75Ly, Lz, sign = -1, A = -A, k,),
        init_vortex_line(; x = 0.75Lx, y = 0.25Ly, Lz, sign = -1, A = -A, k,),
        init_vortex_line(; x = 0.75Lx, y = 0.75Ly, Lz, sign = +1, A = +A, k,),
    ]

    N = 36
    filaments = map(lines) do line
        (; tlims, S, offset,) = line
        dζ = (tlims[2] - tlims[1]) / N
        ζs_base = range(tlims[1] + dζ / 2, tlims[2]; step = dζ)
        @assert length(ζs_base) == N && last(ζs_base) ≈ tlims[2] - dζ / 2
        Filaments.init(ClosedFilament, S.(ζs_base), CubicSplineMethod(); offset)
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
    dt_base = 1.5 * ℓ_min^2 / β        # this criterion seems to provide stability when using RK4
    dt = T_kw / round(T_kw / dt_base)  # make sure that we exactly sample t = T_kw

    fs = copy.(filaments)
    tspan = (0.0, 4.2 * T_kw)

    jprobe = 5
    X_probe = eltype(fs[1])[]
    times = Float64[]
    function callback(iter)
        push!(times, iter.t)
        push!(X_probe, iter.fs[1][jprobe])
    end

    prob = @inferred VortexFilamentProblem(fs, tspan, params_bs)
    # @inferred BiotSavart.init_cache(prob.p; timer = TimerOutput("VF"))

    iter = @inferred init(
        prob, RK4();
        dt, refinement = NoRefinement(),
        callback,
    )
    # solve!(iter)

    # Analyse probed data


    times, X_probe
end
