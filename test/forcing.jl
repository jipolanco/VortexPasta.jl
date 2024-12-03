using VortexPasta
using VortexPasta.Forcing
using StaticArrays
using LinearAlgebra
using Random
using Test

function check_forcing(forcing::FourierNormalFluidForcing)
    (; data, v_rms,) = forcing
    (; cs, qs, Δks,) = data
    N = length(qs[1])  # number of dimensions (usually 3)
    variance = zero(v_rms)
    divergence = zero(v_rms)
    for i ∈ eachindex(qs, cs)
        k⃗ = SVector(qs[i] .* Δks)
        k² = sum(abs2, k⃗)
        @test k² > 0  # make sure zero mode is not included (even if we put kmin = 0.0)
        @test k⃗[1] ≥ 0  # Hermitian symmetry: modes kx < 0 are not included
        u⃗ = cs[i]
        u² = sum(abs2, u⃗)
        factor = iszero(k⃗[1]) ? 1 : 2  # Hermitian symmetry
        divergence += abs2(k⃗ ⋅ u⃗)
        variance += factor * u²
    end
    @test divergence < eps(variance)^2  # divergence is zero
    @test N * v_rms^2 ≈ variance        # check variance
    nothing
end

@testset "Normal fluid forcing" begin
    T = Float32
    rng = Xoshiro(42)
    Ls = T.((2π, 2π, 8π))  # elongated domain
    forcing = @inferred ConstantFourierNormalFluidForcing(rng, Ls; α = 0.2, kmin = 0.0, kmax = 2.5, v_rms = 1.0)
    @test forcing.v_rms isa T
    check_forcing(forcing)
end
