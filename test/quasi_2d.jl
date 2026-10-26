using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.PredefinedCurves
using VortexPasta.VectorsOfVectors: VectorOfVectors
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using LinearAlgebra: norm
using StableRNGs: StableRNG
using Test

function init_params_biot_savart(splitting::AbstractEwaldSplitting)
    ParamsBiotSavart(
        Float64;
        Γ = 2.0, a = 1e-8,
        splitting,
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4))  # ~1e-6 accuracy
    )
end

##

Np = 4
method = CubicSplineMethod()
L = 2π
Ls = (L, L, L)

# Create random array of aligned vortices with identical helical perturbations.
nvort = 10
rng = StableRNG(42)

fs = map(1:nvort) do n
    p = PeriodicLine()
    x = rand(rng) * Ls[1]
    y = rand(rng) * Ls[2]
    orientation = n <= nvort ÷ 2 ? -1 : +1  # neutral charge (not really needed)
    S = define_curve(p; scale = L, translate = (x, y, Ls[3]/2), orientation)
    ζs = range(0, 1; length = 2 * Np + 1)[2:2:end]
    Filaments.init(S, ClosedFilament, ζs, method)
end

Ngrid = 32
splitting = KaiserBesselSplitting(; rtol = 1e-8, Ls, Ns = (Ngrid, Ngrid, Ngrid))
params = init_params_biot_savart(splitting)

tmax = 1.0
prob = VortexFilamentProblem(fs, tmax, params)

dt_kw = BiotSavart.kelvin_wave_period(params, L / Np)
dt = 4 * dt_kw  # this would be unstable without quasi_2D = true

@testset "Scheme: $scheme" for scheme in (RK4(), Strang())
    @testset "fast_term = $fast_term" for fast_term in (LocalTerm(), ShortRangeTerm())  # not really relevant for RK4 (only for splitting schemes)
        iter = init(
            prob, scheme;
            dt,
            quasi_2D = true,
            fast_term,
        )
        E₀ = Diagnostics.kinetic_energy(iter)
        for _ in 1:10
            # @show iter.fs[1][1, CurvatureScalar()]
            # Timestepping.save_checkpoint("quasi2d_$(iter.nstep).vtkhdf", iter) do io
            #     io["velocity"] = iter.vs
            #     io["streamfunction"] = iter.ψs
            #     io["tangent"] = UnitTangent()
            # end
            step!(iter)
        end
        ρ_max = maximum(iter.fs) do f
            maximum(i -> f[i, CurvatureScalar()], eachindex(f))
        end
        # @show ρ_max
        if scheme === RK4()
            @test ρ_max < 1e-15  # good!!
        else
            @test ρ_max < 0.01  # can we do better with splitting scheme??
        end
        E = Diagnostics.kinetic_energy(iter)
        # @show (E - E₀) / E₀
        @test isapprox(E, E₀; rtol = 1.0e-2)  # simulation didn't explode
    end
end
