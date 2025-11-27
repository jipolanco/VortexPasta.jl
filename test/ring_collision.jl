# This tests a simple case with more than one vortex filament.

using Test
using LinearAlgebra: norm, normalize, ⋅
using StaticArrays
using Statistics: mean, std
using OpenCL, pocl_jll
using JET: JET
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Diagnostics: Diagnostics

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

function init_ring_filament(; R, z, sign)
    tlims = (0.0, 2.0)
    S(t) = SVector(π + R * cospi(sign * t), π + R * sinpi(sign * t), z)
    (; R, z, sign, tlims, S,)
end

function test_ring_collision(;
        backend_long = NonuniformFFTsBackend(),
        backend_short = CellListsBackend(),
    )
    R = π / 3  # ring radius
    L = π / 8  # ring distance

    rings = [
        init_ring_filament(; R, z = π - L / 2, sign = +1),
        init_ring_filament(; R, z = π + L / 2, sign = -1),
    ]

    filaments = map(rings) do ring
        (; tlims, S,) = ring
        N = 250  # large enough for autotuning to do some work
        ζs = range(tlims...; length = N + 1)[1:N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    end

    Γ = 2.0
    a = 1e-6
    Δ = 1/4
    Lbox = 3π
    β = 3.5  # accuracy parameter

    params = BiotSavart.autotune(
        filaments, β;
        Γ, a, Δ, Ls = Lbox,
        backend_long,
        backend_short,
        verbose = VERBOSE,
        Cstart = 1.5,
    )
    # println(params)

    cache = BiotSavart.init_cache(params)
    vs = map(f -> similar(nodes(f)), filaments)
    ψs = map(similar, vs)
    fields = (velocity = vs, streamfunction = ψs)
    compute_on_nodes!(fields, cache, filaments)

    spectrum = let
        ks, Ek = Diagnostics.energy_spectrum(cache)
        (; ks, Ek,)
    end

    n = 0
    for (f, v) ∈ zip(filaments, vs)
        velocity_perpendicular_to_tangent = true
        velocity_towards_outside = true
        for i ∈ eachindex(f)
            n += 1
            Ẋ = f[i, Derivative(1)]
            Ẍ = f[i, Derivative(2)]
            if i in (firstindex(f), 4, lastindex(f))  # arbitrary
                @test Ẋ == cache.shortrange.pointdata.derivatives_on_nodes[1][n]
                @test Ẍ == cache.shortrange.pointdata.derivatives_on_nodes[2][n]
            end
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

    JET.@test_opt Diagnostics.stretching_rate(filaments[1], vs[1])
    JET.@test_opt Diagnostics.stretching_rate(filaments[1], vs[1]; quad = GaussLegendre(2))
    JET.@test_call Diagnostics.stretching_rate(filaments[1], vs[1])
    JET.@test_call Diagnostics.stretching_rate(filaments[1], vs[1]; quad = GaussLegendre(2))

    for (f, v) ∈ zip(filaments, vs)
        v_perp = map(v⃗ -> norm(setindex(v⃗, 0.0, 3)), v)
        @test std(v_perp) < 1e-4
        @test mean(v_perp) ≈ vradial_expected rtol=0.002
        dLdt = Diagnostics.stretching_rate(f, v)
        @test dLdt ≈ stretching_rate_expected rtol=0.002
    end

    (; fields..., to = cache.to, spectrum,)
end

@testset "Ring collision" begin
    cpu = test_ring_collision(
        backend_long = NonuniformFFTsBackend(CPU()),
        backend_short = CellListsBackend(CPU()),
    )
    backend = OpenCLBackend()
    gpu = test_ring_collision(
        backend_long = NonuniformFFTsBackend(backend),
        backend_short = CellListsBackend(backend),
    )
    # Note: due to autotuning (which is quite random), the two sets of results can be a bit
    # different, thus we use a relative large `rtol`.
    @testset "Compare CPU / $backend" begin
        @testset "Velocity and streamfunction" begin
            @test cpu.velocity ≈ gpu.velocity rtol=1e-6
            @test cpu.streamfunction ≈ gpu.streamfunction rtol=1e-6
        end
        @testset "Energy spectra" begin
            # Compare energy spectra (first 8 modes only, since kmax may differ between both cases
            # due to autotuning). We make sure we have computed enough modes, but this is basically
            # almost the case.
            if length(cpu.spectrum.ks) ≥ 9 && length(gpu.spectrum.ks) ≥ 9
                @test cpu.spectrum.ks[1:8] ≈ gpu.spectrum.ks[1:8]
                @test cpu.spectrum.Ek[1:8] ≈ gpu.spectrum.Ek[1:8] rtol=1e-5
            end
        end
    end
end
