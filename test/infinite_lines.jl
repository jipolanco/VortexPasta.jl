# This tests infinite (non-closed) filaments with different discretisation methods.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.BiotSavart

# Initialise nearly straight vortex line with helicoidal perturbation.
function init_vortex_line(; x, y, Lz = 2π, sign, A = 0.0, k::Int = 1)
    tlims = (-0.5, 0.5)
    S(t) = SVector(
        x + A * cospi(4π * k * t / Lz),
        y + A * sinpi(4π * k * t / Lz),
        (0.5 + sign * t) * Lz,
    )
    offset = setindex(zero(S(0.0)), sign * Lz, 3)
    (; x, y, Lz, sign, tlims, S, offset,)
end

function test_infinite_lines(method)
    lines = let
        A = 0.08
        k = 2
        h = π / 2
        [
            init_vortex_line(; x = 1h, y = 1h, Lz = 2π, sign = +1, A, k,),
            init_vortex_line(; x = 1h, y = 3h, Lz = 2π, sign = -1, A, k,),
            init_vortex_line(; x = 3h, y = 1h, Lz = 2π, sign = -1, A, k,),
            init_vortex_line(; x = 3h, y = 3h, Lz = 2π, sign = +1, A, k,),
        ]
    end

    N = 32
    filaments = map(lines) do line
        (; tlims, S, offset,) = line
        ζs = range(tlims...; length = N + 1)[1:N]
        @inferred Filaments.init(ClosedFilament, S.(ζs), method; offset)
    end

    @testset "Filaments" begin
        f = first(filaments)
        update_coefficients!(f)
        Xoffset = end_to_end_offset(f)
        @test Xoffset[1] == Xoffset[2] == 0
        @test abs(Xoffset[3]) == 2π
        # Xoffset = Vec3(0, 0, 2π)
        ta, tb = knotlims(f)
        T = tb - ta
        t₀ = 0.1  # arbitrary location
        @test f(t₀) + Xoffset ≈ f(t₀ + T)
        @test f(t₀, Derivative(1)) ≈ f(t₀ + T, Derivative(1))
        @test f(t₀, Derivative(2)) ≈ f(t₀ + T, Derivative(2))
    end

    @testset "Change offset" begin
        f = first(filaments)
        fc = copy(f)
        off = Vec3(1, 1, 1)
        g = @inferred Filaments.change_offset(fc, off)
        for (u, v) ∈ zip(Filaments.allvectors.((fc, g))...)
            @test u === v  # we're reusing the same arrays (for nodes, coefficients, etc...)
        end
        update_coefficients!(g)
        @test end_to_end_offset(g) == off
        ta, tb = knotlims(g)
        T = tb - ta
        t₀ = 0.1  # arbitrary location
        @test g(t₀) + off ≈ g(t₀ + T)
    end

    Ls = (2π, 2π, 2π)
    Ns = (1, 1, 1) .* 64
    kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
    α = kmax / 4
    rcut = 4 / α

    params = ParamsBiotSavart(;
        Γ = 2.0,
        a = 1e-6,
        Δ = 1/4,
        Ls, Ns, rcut, α,
    )
    cache = @inferred BiotSavart.init_cache(params, filaments)
    vs = map(f -> similar(nodes(f)), filaments)
    velocity_on_nodes!(vs, cache, filaments)

    @testset "Velocities" begin
        vnorms = norm.(first(vs))
        vmean = mean(vnorms)
        vstd = std(vnorms)
        # The velocity std should be approximately zero.
        # This condition fails for continuity == 0 (HermiteInterpolation(0), i.e. linear segments).
        continuity = Filaments.continuity(Filaments.interpolation_method(first(filaments)))
        if continuity ≥ 1
            @test vstd / abs(vmean) < 1e-3
        end
    end

    nothing
end

@testset "Infinite lines" begin
    methods = (
        "FiniteDiff(2) / Hermite(2)" => FiniteDiffMethod(2, HermiteInterpolation(2)),
        "FiniteDiff(2) / Hermite(1)" => FiniteDiffMethod(2, HermiteInterpolation(1)),
        "FiniteDiff(2) / Hermite(0)" => FiniteDiffMethod(2, HermiteInterpolation(0)),
        "CubicSpline" => CubicSplineMethod(),
    )
    @testset "$name" for (name, method) ∈ methods
        test_infinite_lines(method)
    end
end
