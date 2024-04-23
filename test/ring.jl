using Test
using LinearAlgebra: norm, normalize, ⋅
using Statistics: mean, std
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Quadratures
using VortexPasta.Diagnostics: Diagnostics
using Random
using StableRNGs: StableRNG

function init_ring_filament(
        N::Int, R = π / 3;
        noise = 0.0, method = CubicSplineMethod(), rng = nothing,
        kwargs...,
    )
    S = define_curve(Ring(); scale = R, translate = π)
    tlims = (0, 1)
    ζs_base = range(tlims...; length = 2N + 1)[2:2:end]
    dζ = step(ζs_base)
    ζs = collect(ζs_base)
    @assert length(ζs) == N
    if !iszero(noise)
        rng_ = rng === nothing ? StableRNG(42) : rng
        ζs .+= noise * dζ * rand(rng_, N)  # `noise` should be ∈ ]-1, 1[
    end
    f = @inferred Filaments.init(S, ClosedFilament, ζs, method; kwargs...)
    (; f, R, noise, method,)
end

# Non-local self-induced velocity of a vortex ring.
# Corresponds to the self-induced velocity of a vortex ring, excluding the
# velocity induced at s⃗ by the local vortex elements [s⃗ - l₋, s⃗ + l₊].
# Here ℓ = √(l₋ l₊) / 4. It is assumed that l₋, l₊ ≪ R.
vortex_ring_nonlocal_velocity(Γ, R, ℓ) = Γ / (4π * R) * log(R / ℓ)

# Total self-induced vortex ring velocity.
vortex_ring_velocity(Γ, R, a; Δ) = Γ / (4π * R) * (log(8R / a) - Δ)

vortex_ring_streamfunction(Γ, R, a; Δ) = Γ / 2π * (log(8R / a) - (1 + Δ))

# Test vortex ring without periodic BCs (i.e. setting α = 0 in Ewald's method, disabling long-range part)
function test_vortex_ring_nonperiodic(ring; quad = GaussLegendre(4))
    (; R, f, noise, method,) = ring
    Γ = 2.4
    a = 1e-6
    Δ = 1/4
    ps = (;
        Γ, a, Δ,
        Ls = Infinity(),
        α = Zero(),
        backend_short = NaiveShortRangeBackend(),
        quadrature = quad,
        lia_segment_fraction = 0.4,
    )
    nquad = length(quad)  # number of quadrature points per segment

    fs = [f]
    params = @inferred ParamsBiotSavart(; ps...)
    cache = @inferred BiotSavart.init_cache(params, fs)
    vs_all = map(similar ∘ nodes, fs)
    ψs_all = map(similar ∘ nodes, fs)
    fields = (velocity = vs_all, streamfunction = ψs_all)

    compute_on_nodes!(fields, cache, fs)
    vs = only(vs_all)  # velocity of the first (and only) filament
    ψs = only(ψs_all)

    @testset "Total velocity" begin
        v⃗_mean = mean(vs)
        U = norm(v⃗_mean)
        # All points have basically the same velocity (since there's no periodicity to introduce anisotropy).
        # This assumes that the point distribution is uniform; otherwise the variance is much higher.
        # @test all(std(vs) .< U * 1e-12)
        U_expected = vortex_ring_velocity(ps.Γ, R, ps.a; Δ = ps.Δ)
        # @show method noise quad f.parametrisation (U - U_expected) / U_expected U U_expected
        rtol = if quad === NoQuadrature()
            0.1  # the error can be huge!
        elseif quad === GaussLegendre(1)
            7e-3
        elseif method isa CubicSplineMethod
            3e-3
        elseif f.parametrisation === CentripetalParametrisation()
            2e-3
        else
            1e-4
        end
        @test isapprox(U, U_expected; rtol)  # the tolerance will mainly depend on the vortex resolution N
    end

    @testset "Total streamfunction" begin
        ψs_para = map(eachindex(ψs)) do i
            t̂ = f[i, UnitTangent()]
            ψs[i] ⋅ t̂
        end
        ψs_perp = map(eachindex(ψs)) do i
            t̂ = f[i, UnitTangent()]
            norm(ψs[i] - ψs_para[i] * t̂)
        end
        # Check that ψ⃗ is always parallel to the tangent vector s′ (⇒ ψs_perp ≈ 0).
        err_perp = maximum(splat((a, b) -> abs(a / b)), zip(ψs_perp, ψs_para))
        tol_noise = noise / 100
        if f.parametrisation === CentripetalParametrisation()
            tol_noise *= 2
        end
        @test err_perp < max(tol_noise, 1e-10)  # use `max` in case noise == 0
        ψ_mean = mean(ψs_para)
        # @test std(ψs_para) < 2 * eps(ψ_mean)  # ψ ⋅ s⃗ along the curve doesn't vary
        ψ_expected = vortex_ring_streamfunction(ps.Γ, R, ps.a; Δ = ps.Δ)
        # @show quad (ψ_mean - ψ_expected) / ψ_expected
        rtol = nquad == 1 ? 0.02 : 1e-4
        @test isapprox(ψ_mean, ψ_expected; rtol)
    end

    # Estimate normalised energy
    @testset "Normalised energy" begin
        Ls = BiotSavart.periods(params)  # all Infinity() in this case
        E_expected = (Γ^2 * R/2) * (log(8R / a) - (Δ + 1))  # energy per unit density
        E_estimated = Diagnostics.kinetic_energy_from_streamfunction(ψs, f, Γ, Ls)  # energy per unit density
        # E_alt = Diagnostics.kinetic_energy_nonperiodic(vs, f, Γ)
        # @show E_expected E_estimated E_alt
        # @show f.parametrisation (E_expected - E_estimated) / E_expected
        rtol = nquad == 1 ? 0.02 : 1e-4
        if f.parametrisation === CentripetalParametrisation()
            rtol = max(rtol, 2e-3)
        end
        @test isapprox(E_expected, E_estimated; rtol)
    end

    # Check velocity excluding LIA (i.e. only non-local integration)
    @testset "Non-local velocity" begin
        i = 1
        integrand(f, i, ζ) = norm(f(i, ζ, Derivative(1)))  # for segment length
        ℓ₋ = integrate(integrand, f, i - 1, GaussLegendre(8))
        ℓ₊ = integrate(integrand, f, i - 0, GaussLegendre(8))
        ℓ = sqrt(ℓ₋ * ℓ₊) / 4 * ps.lia_segment_fraction
        BiotSavart.velocity_on_nodes!([vs], cache, [f]; LIA = Val(false))
        U_nonlocal_expected = vortex_ring_nonlocal_velocity(ps.Γ, R, ℓ)
        U_nonlocal = norm(vs[i])
        rtol = if quad === NoQuadrature()
            0.3  # the error can be huge!!
        elseif quad === GaussLegendre(1)
            0.03
        elseif f.parametrisation === CentripetalParametrisation()
            0.01
        else
            max(1e-3, noise / 400)
        end
        # @show f.parametrisation (U_nonlocal - U_nonlocal_expected) / U_nonlocal_expected rtol
        @test isapprox(U_nonlocal, U_nonlocal_expected; rtol)
    end

    nothing
end

# Check convergence of LIA using quadratures for better accuracy.
function test_local_induced_approximation(ring)
    (; R, f, noise, method,) = ring
    i = firstindex(f)
    ps = (;
        a = 1e-6,
        Δ = 1/4,
        Γ = 4.2,
        segment_fraction = 0.2,
    )
    quad = GaussLegendre(8)  # for accurate estimation of arc length
    lims = @inferred BiotSavart.lia_integration_limits(ps.segment_fraction)
    arclength(j, limits) = integrate((f, j, ζ) -> norm(f(j, ζ, Derivative(1))), f, j, quad; limits)
    ℓ₋ = arclength(i - 1, lims[1])
    ℓ₊ = arclength(i, lims[2])
    ℓ = sqrt(ℓ₋ * ℓ₊) / 4
    v_expected = vortex_ring_velocity(ps.Γ, R, ps.a; ps.Δ) - vortex_ring_nonlocal_velocity(ps.Γ, R, ℓ)
    v_base = norm(
        BiotSavart.local_self_induced_velocity(f, i; quad = nothing, ps...),
    )  # without quadrature
    v_quad = map(1:8) do n
        quad = GaussLegendre(n)
        norm(BiotSavart.local_self_induced_velocity(f, i; quad, ps...))
    end

    # Things converge quite quickly; in this case GaussLegendre(1) seems to be enough.
    # @show method, noise
    # @show (v_base - v_expected) / v_expected
    # @show (v_quad[1] .- v_expected) ./ v_expected
    # @show (v_quad[2] .- v_expected) ./ v_expected
    # @show (v_quad[3] .- v_expected) ./ v_expected

    @test isapprox(v_expected, v_base; rtol = 5e-3)

    rtol = if method isa FourierMethod
        4e-12
    elseif method isa QuinticSplineMethod
        noise == 0 ? 3e-6 : 3e-4
    elseif method isa CubicSplineMethod
        noise == 0 ? 4e-3 : 5e-3
    end
    @test isapprox(v_expected, v_quad[1]; rtol)
    @test isapprox(v_expected, v_quad[2]; rtol)
    @test isapprox(v_expected, v_quad[3]; rtol)

    @testset "Fit circle" begin
        # Alternative estimation by fitting a circle (as in Schwarz PRB 1985).
        # In the specific case of a vortex ring, this should give a perfect
        # estimation of the curvature and binormal vectors.
        v⃗_base_circle = BiotSavart.local_self_induced_velocity(
            f, i; quad = nothing, fit_circle = true, ps...,
        )
        @test normalize(v⃗_base_circle) ≈ Base.setindex(zero(v⃗_base_circle), 1, 3)  # has the right direction
        @test isapprox(norm(v⃗_base_circle), v_expected; rtol = 1e-3)
    end

    nothing
end

function biot_savart_parameters(; L, β = 4, α = 24/L, kws...)
    rcut = β / α
    if L === Infinity()
        M = Zero()
    else
        kmax_wanted = β * 2α
        M = nextprod((2, 3), (L / π) * kmax_wanted)
    end
    ParamsBiotSavart(;
        Γ = 1.0,
        α,
        Ls = (L, L, L),
        Ns = (M, M, M),
        rcut,
        a = 1e-8, Δ = 1/4,
        kws...
    )
end

function test_lia_segment_fraction()
    N = 64 * 1
    method = QuinticSplineMethod()
    ring = init_ring_filament(N; method)
    fs = [ring.f]

    L = Infinity()
    params = biot_savart_parameters(; L, lia_segment_fraction = 0.2)
    # println(params)
    cache = BiotSavart.init_cache(params, fs)

    vs_all = map(similar ∘ nodes, fs)
    ψs_all = map(similar ∘ nodes, fs)
    fields = (velocity = vs_all, streamfunction = ψs_all)
    compute_on_nodes!(fields, cache, fs)

    # seg = Filaments.Segment(fs[1], 3)
    # x⃗ = fs[1](3, 0.3)
    # @time y⃗ = BiotSavart.integrate_biot_savart(
    #     BiotSavart.Velocity(),
    #     BiotSavart.FullIntegrand(),
    #     seg, x⃗,
    #     params.common,
    # )
    # @time y⃗ = BiotSavart.integrate_biot_savart(
    #     BiotSavart.Velocity(),
    #     BiotSavart.FullIntegrand(),
    #     seg, x⃗,
    #     params.common,
    # )

    vz = map(v -> v.z, vs_all[1])
    vz_avg = mean(vz)
    vz_std = std(vz; corrected = false)
    # @show vz_avg vz_std/vz_avg
    @test vz_std/vz_avg < 1e-12

    # This is valid in the non-periodic case only
    v_expected = vortex_ring_velocity(params.Γ, ring.R, params.a; params.Δ)
    err_v = abs(vz_avg - v_expected) / v_expected
    # @show err_v
    @test err_v < 2e-6  # requires using lia_segment_fraction

    nothing
end

@testset "Vortex ring velocity" begin
    N = 32  # number of discretisation points
    # `noise`: non-uniformity of filament discretisation (noise == 0 means uniform discretisation)
    for method ∈ (CubicSplineMethod(), QuinticSplineMethod())
        @testset "$method, N = $N, noise = $noise" for noise ∈ (0.0, 0.2)
            ring = init_ring_filament(N; noise, method)
            @testset "Non-periodic" begin
                # Note: GaussLegendre(1) and NoQuadrature() should give nearly the same results.
                quads = (GaussLegendre(4), GaussLegendre(1), NoQuadrature())
                @testset "Quadrature: $quad" for quad ∈ quads
                    test_vortex_ring_nonperiodic(ring; quad)
                end
            end
            @testset "LIA" begin
                test_local_induced_approximation(ring)
            end
        end
    end
    @testset "CentripetalParametrisation" begin
        local parametrisation = CentripetalParametrisation()
        local method = QuinticSplineMethod()
        local noise = 0.2
        local ring = init_ring_filament(N; noise, method, parametrisation)
        test_vortex_ring_nonperiodic(ring; quad = GaussLegendre(4))
    end
    method = FourierMethod()
    @testset "$method" begin
        ring = init_ring_filament(N; noise = 0.0, method)
        @testset "Non-periodic" begin
            test_vortex_ring_nonperiodic(ring; quad = GaussLegendre(1))
        end
        @testset "LIA" begin
            test_local_induced_approximation(ring)
        end
    end
    @testset "LIA segment fraction" begin
        test_lia_segment_fraction()
    end
end
