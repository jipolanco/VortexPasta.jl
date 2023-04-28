using Test
using LinearAlgebra: norm
using Statistics: mean, std
using VortexPasta.Filaments
using VortexPasta.BiotSavart

function init_ring_filament(N::Int, R = π / 3)
    S(t) = π .+ R .* Vec3(cospi(t), sinpi(t), zero(t))
    tlims = (0, 2)
    ζs = collect(range(tlims...; length = N + 1)[1:N])
    # ζs .+= rand(N) / 2N
    f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    (; f, R,)
end

# Non-local self-induced velocity of a vortex ring.
# Corresponds to the self-induced velocity of a vortex ring, excluding the
# velocity induced at s⃗ by the local vortex elements [s⃗ - l₋, s⃗ + l₊].
# Here ℓ = √(l₋ l₊) / 4. It is assumed that l₋, l₊ ≪ R.
vortex_ring_nonlocal_velocity(Γ, R, ℓ) = Γ / (4π * R) * log(R / ℓ)

# Total self-induced vortex ring velocity.
vortex_ring_velocity(Γ, R, a; Δ) = Γ / (4π * R) * (log(8R / a) - Δ)

# Test vortex ring without periodic BCs (i.e. setting α = 0 in Ewald's method, disabling long-range part)
function test_vortex_ring_nonperiodic(ring)
    (; R, f,) = ring
    N = length(f)
    ps = (;
        Γ = 2.4,
        a = 1e-6,
        Δ = 1/4,
        Ls = (∞, ∞, ∞),
        α = Zero(),
        quadrature_short = GaussLegendreQuadrature(4),
    )
    params = @inferred ParamsBiotSavart(; ps...)
    cache = @inferred BiotSavart.init_cache(params)
    vs = similar(nodes(f))
    @testset "Total velocity" begin
        velocity_on_nodes!(vs, cache, f)
        @show vs[1]
        v⃗_mean = mean(vs)
        U = norm(v⃗_mean)
        # All points have basically the same velocity (since there's no periodicity to introduce anisotropy).
        # This assumes that the point distribution is uniform; otherwise the variance is much higher.
        @test all(std(vs) .< U * 1e-12)
        U_expected = vortex_ring_velocity(ps.Γ, R, ps.a; Δ = ps.Δ)
        @show (U - U_expected) / U_expected
        @test isapprox(U, U_expected; rtol = 1e-4)  # the tolerance will mainly depend on the vortex resolution N
    end
    # Check velocity excluding LIA (i.e. only non-local integration)
    @testset "Non-local velocity" begin
        ts = knots(f)
        i = 1
        ℓ = sqrt((ts[i + 1] - ts[i]) * (ts[i] - ts[i - 1])) / 4  # = √(l₋ l₊) / 4
        fill!(vs, zero(eltype(vs)))
        BiotSavart.add_short_range_velocity_self!(vs, cache.shortrange, f; LIA = false)
        U_nonlocal_expected = vortex_ring_nonlocal_velocity(ps.Γ, R, ℓ)
        U_nonlocal = norm(vs[i])
        @show (U_nonlocal - U_nonlocal_expected) / U_nonlocal_expected
        @test isapprox(U_nonlocal, U_nonlocal_expected; rtol = 1e-3)
    end
    nothing
end

# Check convergence of LIA using quadratures for better accuracy.
function test_local_induced_approximation(ring)
    (; R, f,) = ring
    i = firstindex(f)
    ps = (;
        a = 1e-6,
        Δ = 1/4,
        Γ = 4.2,
    )
    quad = GaussLegendreQuadrature(8)  # for accurate estimation of arc length
    arclength(j) = integrate(ζ -> norm(f(j, ζ, Derivative(1))), f, j, quad)
    ℓ₋ = arclength(i - 1)
    ℓ₊ = arclength(i)
    ℓ = sqrt(ℓ₋ * ℓ₊) / 4
    v_expected = vortex_ring_velocity(ps.Γ, R, ps.a; ps.Δ) - vortex_ring_nonlocal_velocity(ps.Γ, R, ℓ)
    v_base = norm(BiotSavart.local_self_induced_velocity(f, i; quad = nothing, ps...))  # without quadrature
    v_quad = map(1:8) do n
        quad = GaussLegendreQuadrature(n)
        norm(BiotSavart.local_self_induced_velocity(f, i; quad, ps...))
    end
    # Things converge quite quickly; in this case GaussLegendreQuadrature(2) seems to be enough.
    @show (v_base - v_expected) / v_expected
    @show (v_quad .- v_expected) ./ v_expected
    @test isapprox(v_expected, v_base; rtol = 1e-2)
    @test isapprox(v_expected, v_quad[1]; rtol = 1e-2)
    @test isapprox(v_expected, v_quad[2]; rtol = 6e-6)
    @test isapprox(v_expected, v_quad[3]; rtol = 3e-6)
    nothing
end

@testset "Vortex ring" begin
    N = 32
    ring = init_ring_filament(N)
    @testset "Non-periodic" begin
        test_vortex_ring_nonperiodic(ring)
    end
    @testset "LIA" begin
        test_local_induced_approximation(ring)
    end
end
