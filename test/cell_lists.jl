using Test
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using Rotations: AngleAxis
using LinearAlgebra: I
using StaticArrays

function vortex_ring(; R, origin = π * Vec3(1, 1, 1), sign, rotation = I,)
    tlims = (0.0, 2.0)
    S(t) = origin + rotation * SVector(R * cospi(sign * t), R * sinpi(sign * t), 0)
    (; R, origin, sign, tlims, S, offset = zero(Vec3{Float64}),)
end

@testset "SegmentCellList" begin
    Ls = (2π, 2π, 2π)
    ncells_expected = 8
    rcut = 1.01 * 2π / ncells_expected

    # Linked vortex rings
    rings = [
        vortex_ring(; R = 0.6π, origin = Vec3(π - 0.3π, π, π), sign = +1),
        vortex_ring(; R = 0.6π, origin = Vec3(π + 0.3π, π, π), rotation = AngleAxis(π/2, 1, 0, 0), sign = -1),
    ]
    N = 32

    fs = map(rings) do ring
        (; tlims, S, offset,) = ring
        dζ = (tlims[2] - tlims[1]) / N
        ζs = range(tlims[1] + dζ / 2, tlims[2]; step = dζ)
        @assert length(ζs) == N && last(ζs) ≈ tlims[2] - dζ / 2
        @assert step(ζs) ≈ dζ
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod(); offset)
    end

    cl = @inferred BiotSavart.SegmentCellList(fs, rcut, Ls)
    @test size(cl.segments) == (ncells_expected, ncells_expected, ncells_expected)
    @test sum(length, cl.segments) == sum(f -> length(segments(f)), fs)  # contains all segments
end
