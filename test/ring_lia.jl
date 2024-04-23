using Test
using VortexPasta.Filaments
using VortexPasta.PredefinedCurves
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using UnicodePlots: UnicodePlots, lineplot, lineplot!
using Statistics: mean, std
using Random: Random
using StableRNGs: StableRNG
using Rotations: Rotations

function generate_biot_savart_parameters(::Type{T}; periodic = Val(true)) where {T}
    # Parameters are relevant for HeII if we interpret dimensions in cm and s.
    Γ = 9.97e-4
    a = 1e-8
    Δ = 1/4
    if periodic === Val(true)
        L = 2π
        Ngrid = 32
        backend_short = CellListsBackend(2)
    else
        L = Infinity()
        Ngrid = 0
        backend_short = NaiveShortRangeBackend()
    end
    Ls = (L, L, L)
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    α = kmax / 6
    rcut = 5 / α
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short,
        backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4)),
        quadrature = GaussLegendre(2),
        lia_segment_fraction = 0.1,
    )
end

function generate_ring(
        ::Type{T}, N, method = QuinticSplineMethod();
        R, noise = 0.0,
    ) where {T}
    p = Ring()
    scale = R
    rng = StableRNG(4242)
    ts = collect(range(T(0), T(1); length = N + 1)[1:N])
    for i ∈ eachindex(ts)
        ts[i] += rand(rng, T) * noise / N
    end
    S = define_curve(p; scale,)
    [
        Filaments.init(S, ClosedFilament{T}, ts, method),
    ]
end

function test_lia_ring(
        ::Type{T}, N, method;
        scheme = RK4(),
        R, periodic = Val(true),
        noise = 0.0,
        verbose = false,
    ) where {T}
    fs = @inferred generate_ring(T, N, method; R, noise)
    params = @inferred generate_biot_savart_parameters(T; periodic)
    (; Γ,) = params

    τ = R^2 / Γ
    tmax = 0.98 * τ
    tspan = (0.0, tmax)
    prob = VortexFilamentProblem(fs, tspan, params)

    δ = minimum_node_distance(fs)
    dt = BiotSavart.kelvin_wave_period(params, δ) * 1.2
    iter = init(prob, scheme; dt, fold_periodic = false, LIA = true)

    # Total self-induced vortex ring velocity.
    v_ring = params.Γ / (4π * R) * (log(8R / params.a) - params.Δ)

    vz = getindex.(iter.vs[1], 3)
    vz_mean = mean(vz)
    vz_std = std(vz)
    verbose && @show vz_mean/v_ring vz_std/vz_mean
    @test 0.1 * v_ring < vz_mean < 0.9 * v_ring  # since we only include the local part, vz_mean < v_ring
    if noise == 0
        @test vz_std / vz_mean < eps(T) * 1e3
    else
        @test vz_std / vz_mean < noise * 0.02
    end

    E = @inferred Diagnostics.kinetic_energy(iter; quad = GaussLegendre(2))

    times = [iter.t]
    energy = [E]

    while iter.t < tmax
        step!(iter)
        local (; nstep, t, fs,) = iter
        local f = fs[1]
        if verbose && nstep % 10 == 0
            zs = getindex.(f, 3)
            zmean = mean(zs)
            zstd = std(zs)
            @show nstep, t/τ, zmean/R, zstd/R
        end
        E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(2))
        push!(times, t)
        push!(energy, E)
        # write_vtkhdf("links_$nstep.vtkhdf", fs) do io
        #     io["CurvatureVector"] = CurvatureVector()
        # end
    end

    times_norm = times ./ τ
    Emean = mean(energy)
    Estd = std(energy)
    verbose && @show Estd / Emean
    if noise == 0
        @test Estd / Emean < eps(T) * 50
    else
        @test Estd / Emean < noise * 0.02
    end

    if verbose
        let plt = lineplot(
                times_norm, energy ./ energy[begin];
                xlabel = "t Γ / R²", ylabel = "E / E₀",
                title = "LIA ring",
            )
            println(plt)
        end
        println(iter.to)
    end

    nothing
end

@testset "LIA ring" begin
    N = 48
    method = QuinticSplineMethod()
    R = 1.2
    periodic = Val(false)
    verbose = false
    @testset "T = $T" for T ∈ (Float32, Float64)
        @testset "Noise = $noise" for noise ∈ (0.0, 0.5)
            test_lia_ring(T, N, method; R, periodic, noise, verbose)
        end
    end
    # It doesn't make much sense to use splitting schemes with LIA-only, but we check that
    # they work correctly anyways (in the sense that they set the non-local contributions to
    # 0).
    @testset "Splitting scheme" begin
        local scheme = Strang(RK4())
        test_lia_ring(Float32, N, method; R, periodic, verbose, scheme)
    end
end
