using Test
using VortexPasta.Filaments
using VortexPasta.PredefinedCurves
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using UnicodePlots: UnicodePlots, lineplot, lineplot!
using Rotations: Rotations

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

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
    β = 3.5
    α = kmax / 2β
    rcut = β / α
    backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4))
    # Test unexported functions
    @test BiotSavart.oversampling_factor(backend_long) == 1.5
    @test BiotSavart.half_support(backend_long) == 4
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short,
        backend_long,
        quadrature = GaussLegendre(2),
    )
end

function generate_linked_rings(
        ::Type{T}, N, method = QuinticSplineMethod();
        R,
    ) where {T}
    p = Ring()
    δ = R/2
    scale = R
    rotate = Rotations.RotX(π/2)
    Sl = define_curve(p; scale, translate = (-δ, 0, 0))           # left ring
    Sr = define_curve(p; scale, rotate, translate = (+δ, 0, 0))  # right ring
    [
        Filaments.init(Sl, ClosedFilament{T}, N, method),
        Filaments.init(Sr, ClosedFilament{T}, N, method),
    ]
end

function test_linked_rings(
        ::Type{T}, N, method;
        R, periodic = Val(true),
    ) where {T}
    fs = @inferred generate_linked_rings(T, N, method; R)
    params = @inferred generate_biot_savart_parameters(T; periodic)
    (; Γ,) = params

    τ = R^2 / Γ
    tmax = 0.8 * τ  # reconnection seems to happen exactly around t = τ / 2
    tspan = (0.0, tmax)
    prob = VortexFilamentProblem(fs, tspan, params)

    δ = minimum_node_distance(fs)
    # d_crit = R / 20
    d_crit = δ / 2
    dt = BiotSavart.kelvin_wave_period(params, δ)
    reconnect = ReconnectBasedOnDistance(d_crit)
    # @show R d_crit δ tmax dt
    iter = init(prob, RK4(); dt, reconnect, fold_periodic = false)

    E = @inferred Diagnostics.kinetic_energy(iter; quad = GaussLegendre(2))
    H_no_quad = @inferred Diagnostics.helicity(iter; quad = nothing)
    H = @inferred Diagnostics.helicity(iter; quad = GaussLegendre(2))
    @test E isa T
    @test H_no_quad isa T
    @test H isa T

    if periodic === Val(true)
        # Note: in this case the helicity spectrum decays very slowly, so the sum gives a
        # very rough approximation to the actual helicity.
        ks, Hk = @inferred Diagnostics.helicity_spectrum(iter)
        dk = ks[2] - ks[1]
        println(Hk ./ params.Γ^2)
        H_from_spectrum = sum(Hk) * dk
        # @show H H_from_spectrum
        @test 2 * H_from_spectrum < H < H_from_spectrum < 0
    end

    # @show (H - H_no_quad) / H
    @test isapprox(H, H_no_quad)  # there is really no difference between the two

    # In the case of two unknotted linked rings, the helicity is H = 2 Γ² L where L is the
    # linking number. For two linked rings, |L| = 1.
    linking = H / (2 * Γ^2)

    if periodic === Val(true)
        # When periodicity is enabled, the accuracy seems to be mainly controlled by the
        # splitting parameter α/kmax and by the accuracy of the NUFFTs.
        rtol = 4e-6
    else
        rtol = 2e-9  # without periodicity things are really accurate
    end
    # @show linking + 1
    @test isapprox(linking, -1; rtol)

    times = [iter.t]
    energy = [E]
    helicity = [H]
    t_reconnect = -1.0

    while iter.t < tmax
        step!(iter)
        local (; nstep, t, fs,) = iter
        if nstep % 10 == 0
            # @show nstep, t/τ, length(fs)
        end
        E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(2))
        H = Diagnostics.helicity(iter; quad = GaussLegendre(2))
        push!(times, t)
        push!(energy, E)
        push!(helicity, H)
        if t_reconnect < 0 && length(fs) == 1
            t_reconnect = t
        end
        # write_vtkhdf("links_$nstep.vtkhdf", fs) do io
        #     io["CurvatureVector"] = CurvatureVector()
        # end
    end

    times_norm = times ./ τ

    if VERBOSE
        let plt = lineplot(
                times_norm, energy ./ energy[begin];
                xlabel = "t Γ / R²", ylabel = "E / E₀",
                title = "Hopf link reconnection",
            )
            UnicodePlots.vline!(plt, 0.5)
            println(plt)
        end

        let plt = lineplot(
                times_norm, abs.(helicity) ./ (2 * Γ^2);
                xlabel = "t Γ / R²", ylabel = "|H| / 2Γ²",
                ylim = (0.5, 1.0),
            )
            UnicodePlots.vline!(plt, 0.5)
            println(plt)
        end

        println(iter.to)
    end

    @test 0.49 < t_reconnect / τ < 0.51

    nothing
end

@testset "Linked rings" begin
    T = Float32
    N = 48
    method = QuinticSplineMethod()
    R = 1.2
    @testset "Periodic: $periodic" for periodic in (Val(false), Val(true))
        test_linked_rings(T, N, method; R, periodic)
    end
end
