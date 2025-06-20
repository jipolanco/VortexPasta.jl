using Test
using VortexPasta.Filaments
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot
using VortexPasta.BiotSavart
using VortexPasta.Timestepping

##

function biot_savart_parameters()
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

    ParamsBiotSavart(;
        Γ, a, Δ,
        α, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = NonuniformFFTsBackend(),
        quadrature = GaussLegendre(3),
    )
end

function check_tableaus(scheme::Timestepping.ImplicitExplicitScheme)
    E = Timestepping.tableau_explicit(scheme)
    I = Timestepping.tableau_implicit(scheme)
    # These vectors must be the same in each tableau
    @test E.bs == I.bs
    @test E.cs == I.cs
    for i ∈ axes(E.A, 1)
        # Sum of each row must be consistent with `cs`.
        @test sum(E.A[i, :]) ≈ E.cs[i]
        @test sum(I.A[i, :]) ≈ E.cs[i]
    end
    nothing
end

##

@testset "Butcher tableaus" begin
    schemes = (KenCarp3(), KenCarp4(),)
    @testset "$scheme" for scheme ∈ schemes
        check_tableaus(scheme)
    end
end

@testset "Local/non-local splitting" begin
    ##
    params_bs = @inferred biot_savart_parameters()

    f_init = Filaments.init(
        define_curve(TrefoilKnot()), ClosedFilament, 42, CubicSplineMethod(),
    )
    fs_init = [f_init]
    tspan = (0.0, 1.0)

    prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs);

    scheme = IMEXEuler()

    external_velocity(x⃗, t) = Vec3(0.1, 0.0, 0.0)  # add constant external velocity

    @testset "fast_term = $fast_term" for fast_term ∈ (LocalTerm(), ShortRangeTerm())
        iter = @inferred init(prob, scheme; dt = 0.025, fast_term, external_velocity)
        iter = @inferred init(prob, scheme; dt = 0.025, fast_term, alias_u0 = false, external_velocity)

        (; fs, vs, rhs!,) = iter

        # This is just to make sure that eachindex(fs[i], vs[i]) returns linear indices.
        @testset "Ensure linear indexing" begin
            (; fs, vs,) = iter
            @test eachindex(fs) isa AbstractUnitRange
            @test eachindex(vs) isa AbstractUnitRange
            @test eachindex(fs, vs) isa AbstractUnitRange
            let f = first(fs), v = first(vs)
                @test eachindex(f, v) isa AbstractUnitRange
            end
        end

        # Check that slow component (non-local) + fast component (LIA) == full velocity
        @testset "Velocity splitting" begin
            vs_full = similar(vs)
            vs_slow = similar(vs)
            vs_fast = similar(vs)

            rhs!(vs_full, fs, iter.time.t, iter; component = Val(:full))
            rhs!(vs_slow, fs, iter.time.t, iter; component = Val(:slow))
            rhs!(vs_fast, fs, iter.time.t, iter; component = Val(:fast))

            @test vs ≈ vs_full  # vs is computed when the solver is initialised
            @test vs_slow + vs_fast ≈ vs_full
        end
    end

    ##
end
