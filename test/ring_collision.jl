# This tests a simple case with more than one vortex filament.

using Test
using LinearAlgebra: norm, normalize, ⋅
using StaticArrays
using Statistics: mean, std
using JET: JET
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Diagnostics: Diagnostics

function init_ring_filament(; R, z, sign)
    tlims = (0.0, 2.0)
    S(t) = SVector(π + R * cospi(sign * t), π + R * sinpi(sign * t), z)
    (; R, z, sign, tlims, S,)
end

function test_ring_collision()
    R = π / 3  # ring radius
    L = π / 8  # ring distance

    rings = [
        init_ring_filament(; R, z = π - L / 2, sign = +1),
        init_ring_filament(; R, z = π + L / 2, sign = -1),
    ]

    filaments = map(rings) do ring
        (; tlims, S,) = ring
        N = 64
        ζs = range(tlims...; length = N + 1)[1:N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    end

    Γ = 2.0
    a = 1e-6
    Δ = 1/4
    Lbox = 2π
    Ls = (Lbox, Lbox, Lbox)
    rcut = Lbox / 3
    β = 3.5  # accuracy parameter
    α = β / rcut
    kmax = 2α * β
    Ngrid = ceil(Int, kmax * Lbox / π)
    Ns = (Ngrid, Ngrid, Ngrid)

    params = ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
    )
    # println(params)

    cache = BiotSavart.init_cache(params, filaments)
    vs = map(f -> similar(nodes(f)), filaments)
    velocity_on_nodes!(vs, cache, filaments)

    for (f, v) ∈ zip(filaments, vs)
        velocity_perpendicular_to_tangent = true
        velocity_towards_outside = true
        for i ∈ eachindex(f)
            Ẋ = f[i, Derivative(1)]
            Ẍ = f[i, Derivative(2)]
            t̂, ρ⃗ = normalise_derivatives(Ẋ, Ẍ)  # tangent and curvature vectors
            v⃗ = v[i]
            v⃗_perp = setindex(v⃗, 0.0, 3)  # set vz = 0
            velocity_perpendicular_to_tangent &= abs(v⃗ ⋅ t̂) / norm(v⃗) < 1e-4
            velocity_towards_outside &= isapprox(v⃗_perp ⋅ ρ⃗, -norm(v⃗_perp) * norm(ρ⃗))
        end
        @test velocity_perpendicular_to_tangent
        @test velocity_towards_outside
    end

    # Expected radial velocity of a vortex ring (induced only by the other ring).
    # Note that, since we're periodic, there are some tiny deviations from the
    # expected velocity (< 1% in this case) due to periodic copies of the vortices.
    vradial_expected = let N = 100
        ℓ = L / R
        θs = range(0, 2π; length = N + 1)[1:N]
        dθ = step(θs)
        (Γ * ℓ / (4π * R)) * dθ * sum(θs) do θ
            c = cos(θ)
            r² = 2 * (1 - c) + ℓ^2
            c / sqrt(r²)^3
        end
    end

    # Increase of vortex length dL/dt, can be easily computed for a vortex ring with uniform
    # radial velocity.
    stretching_rate_expected = 2π * vradial_expected

    disable_jet = get(ENV, "JULIA_DISABLE_JET_KA_TESTS", "false") ∈ ("true", "1")  # disable JET tests involving KA kernels
    if !disable_jet
        JET.@test_opt Diagnostics.stretching_rate(filaments[1], vs[1])
        JET.@test_opt Diagnostics.stretching_rate(filaments[1], vs[1]; quad = GaussLegendre(2))
        JET.@test_call Diagnostics.stretching_rate(filaments[1], vs[1])
        JET.@test_call Diagnostics.stretching_rate(filaments[1], vs[1]; quad = GaussLegendre(2))
    end

    for (f, v) ∈ zip(filaments, vs)
        v_perp = map(v⃗ -> norm(setindex(v⃗, 0.0, 3)), v)
        @test std(v_perp) < 1e-4
        @test isapprox(mean(v_perp), vradial_expected; rtol = 0.002)
        dLdt = Diagnostics.stretching_rate(f, v)
        @test isapprox(dLdt, stretching_rate_expected; rtol = 0.002)
    end

    nothing
end

@testset "Ring collision" begin
    test_ring_collision()
end
