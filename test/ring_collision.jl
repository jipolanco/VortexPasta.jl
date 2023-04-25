# This tests a simple case with more than one vortex filament.

using Test
using LinearAlgebra: norm, normalize, ⋅
using StaticArrays
using VortexPasta.Filaments
using VortexPasta.BiotSavart

function init_ring_filament(; R, z, sign)
    tlims = (0.0, 2.0)
    S(t) = SVector(π + R * cospi(sign * t), π + R * sinpi(sign * t), z)
    (; R, z, sign, tlims, S,)
end

function test_ring_collision()
    rings = [
        init_ring_filament(; R = π / 3, z = π - π / 16, sign = +1),
        init_ring_filament(; R = π / 3, z = π + π / 16, sign = -1),
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

    # Each vortex has a velocity which can be decomposed as v = v⟂ + vz, where
    # vz is its own self-induced velocity, and v⟂ is the velocity induced by
    # the other vortex. In this case, v⟂ is orthogonal to vz. It is also
    # orthogonal to the filament tangent and is opposite to the curvature vector,
    # i.e. it is directed towards the exterior of the vortex (such that the
    # vortex tends to grow).

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

    nothing
end

@testset "Ring collision" begin
    test_ring_collision()
end
