using Test
using LinearAlgebra: norm, normalize, ⋅
using Statistics: mean, std
using VortexPasta.PredefinedCurves: define_curve, Ring
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Diagnostics: Diagnostics
using Random

function init_ring_filament(N::Int, R = π / 3; noise = 0.0, rng = nothing)
    S = define_curve(Ring(); scale = R, translate = π)
    tlims = (0, 1)
    ζs_base = range(tlims...; length = 2N + 1)[2:2:end]
    dζ = step(ζs_base)
    ζs = collect(ζs_base)
    @assert length(ζs) == N
    if !iszero(noise)
        rng_ = rng === nothing ? MersenneTwister(42) : rng
        ζs .+= noise * dζ * rand(rng_, N)  # `noise` should be ∈ ]-1, 1[
    end
    f = @inferred Filaments.init(S, ClosedFilament, ζs, CubicSplineMethod())
    (; f, R, noise,)
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
    (; R, f, noise,) = ring
    Γ = 2.4
    a = 1e-6
    Δ = 1/4
    ps = (;
        Γ, a, Δ,
        Ls = Infinity(),
        α = Zero(),
        backend_short = NaiveShortRangeBackend(),
        quadrature_short = quad,
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

    @testset "NaiveSegmentIterator" begin
        it = @inferred BiotSavart.nearby_segments(cache.shortrange, f[2])
        @test it isa BiotSavart.NaiveSegmentIterator
        @test eltype(it) === typeof(first(segments(f))) === typeof(first(it))
        segs = @inferred collect(it)
        @test eltype(it) === eltype(segs)
    end

    @testset "Total velocity" begin
        v⃗_mean = mean(vs)
        U = norm(v⃗_mean)
        # All points have basically the same velocity (since there's no periodicity to introduce anisotropy).
        # This assumes that the point distribution is uniform; otherwise the variance is much higher.
        # @test all(std(vs) .< U * 1e-12)
        U_expected = vortex_ring_velocity(ps.Γ, R, ps.a; Δ = ps.Δ)
        @show (U - U_expected) / U_expected
        rtol = nquad == 1 ? 5e-3 : 1e-4
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
        @test all(zip(ψs_para, ψs_perp)) do (a, b)
            b / a < max(noise / 100, 1e-10)  # use `max` in case noise == 0
        end
        ψ_mean = mean(ψs_para)
        # @test std(ψs_para) < 2 * eps(ψ_mean)  # ψ ⋅ s⃗ along the curve doesn't vary
        ψ_expected = vortex_ring_streamfunction(ps.Γ, R, ps.a; Δ = ps.Δ)
        @show (ψ_mean - ψ_expected) / ψ_expected
        rtol = 4e-3 / nquad
        @test isapprox(ψ_mean, ψ_expected; rtol)
    end

    # Estimate normalised energy
    @testset "Normalised energy" begin
        Ls = 1  # sets domain volume to 1 for energy estimation
        E_expected = (Γ^2 * R/2) * (log(8R / a) - (Δ + 1))
        E_estimated = Diagnostics.kinetic_energy_from_streamfunction(ψs, f, Γ, Ls)
        @show (E_expected - E_estimated) / E_expected
        rtol = 4e-3 / nquad
        @test isapprox(E_expected, E_estimated; rtol)
    end

    # Check velocity excluding LIA (i.e. only non-local integration)
    @testset "Non-local velocity" begin
        ts = knots(f)
        i = 1
        ℓ = sqrt((ts[i + 1] - ts[i]) * (ts[i] - ts[i - 1])) / 4  # = √(l₋ l₊) / 4
        fill!(vs, zero(eltype(vs)))
        BiotSavart.set_filaments!(cache.shortrange, [f])
        BiotSavart.add_short_range_velocity!(vs, cache.shortrange, f; LIA = Val(false))
        U_nonlocal_expected = vortex_ring_nonlocal_velocity(ps.Γ, R, ℓ)
        U_nonlocal = norm(vs[i])
        @show (U_nonlocal - U_nonlocal_expected) / U_nonlocal_expected
        rtol = max(16e-3 / nquad^2, noise / 400)
        @test isapprox(U_nonlocal, U_nonlocal_expected; rtol)
    end

    nothing
end

# Check convergence of LIA using quadratures for better accuracy.
function test_local_induced_approximation(ring)
    (; R, f, noise,) = ring
    i = firstindex(f)
    ps = (;
        a = 1e-6,
        Δ = 1/4,
        Γ = 4.2,
    )
    quad = GaussLegendre(8)  # for accurate estimation of arc length
    arclength(j) = integrate((f, j, ζ) -> norm(f(j, ζ, Derivative(1))), f, j, quad)
    ℓ₋ = arclength(i - 1)
    ℓ₊ = arclength(i)
    ℓ = sqrt(ℓ₋ * ℓ₊) / 4
    v_expected = vortex_ring_velocity(ps.Γ, R, ps.a; ps.Δ) - vortex_ring_nonlocal_velocity(ps.Γ, R, ℓ)
    v_base = norm(BiotSavart.local_self_induced_velocity(f, i; quad = nothing, ps...))  # without quadrature
    v_quad = map(1:8) do n
        quad = GaussLegendre(n)
        norm(BiotSavart.local_self_induced_velocity(f, i; quad, ps...))
    end
    # Things converge quite quickly; in this case GaussLegendre(2) seems to be enough.
    @show (v_base - v_expected) / v_expected
    @show (v_quad[1] .- v_expected) ./ v_expected
    @show (v_quad[2] .- v_expected) ./ v_expected
    @show (v_quad[3] .- v_expected) ./ v_expected
    @test isapprox(v_expected, v_base; rtol = 5e-3)
    @test isapprox(v_expected, v_quad[1]; rtol = 1e-2)
    rtol = max(6e-6, noise * 1e-3)
    @test isapprox(v_expected, v_quad[2]; rtol)
    @test isapprox(v_expected, v_quad[3]; rtol)
    @testset "Fit circle" begin
        # Alternative estimation by fitting a circle (as in Schwarz PRB 1985).
        # In the specific case of a vortex ring, this should give a perfect
        # estimation of the curvature and binormal vectors.
        v⃗_base_circle = BiotSavart.local_self_induced_velocity(f, i; quad = nothing, fit_circle = true, ps...)
        @test normalize(v⃗_base_circle) ≈ Base.setindex(zero(v⃗_base_circle), 1, 3)  # has the right direction
        @test isapprox(norm(v⃗_base_circle), v_expected; rtol = 1e-3)
    end
    nothing
end

N = 32  # number of discretisation points

# `noise`: non-uniformity of filament discretisation (noise == 0 means uniform discretisation)

@testset "Vortex ring (N = $N, noise = $noise)" for noise ∈ (0.0, 0.2)
    ring = init_ring_filament(N; noise)
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

if @isdefined(Makie)
    let fig = Figure()
        ax = Axis3(fig[1, 1])
        filamentplot!(ax, ring.f)
        fig
    end
end
