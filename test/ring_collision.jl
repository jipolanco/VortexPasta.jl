# This tests a simple case with more than one vortex filament.

using Test
using LinearAlgebra: norm, normalize, ⋅
using StaticArrays
using Statistics: mean, std
using VortexPasta.Filaments
using VortexPasta.BiotSavart

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
    Ls = (2π, 2π, 2π)
    Ns = (64, 64, 64)
    kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
    α = kmax / 6
    rcut = 4 / α

    params = ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
    )

    cache = BiotSavart.init_cache(params)
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

    for v ∈ vs
        v_perp = map(v⃗ -> norm(setindex(v⃗, 0.0, 3)), v)
        @test std(v_perp) < 1e-4
        @test isapprox(mean(v_perp), vradial_expected; rtol = 0.002)
    end

    nothing
end

@testset "Ring collision" begin
    test_ring_collision()
end
