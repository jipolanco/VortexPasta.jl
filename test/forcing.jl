using VortexPasta
using VortexPasta.Forcing
using StaticArrays
using LinearAlgebra
using Random
using Test

function check_forcing(forcing::FourierNormalFluidForcing)
    (; data, vn_rms,) = forcing
    (; cs, ks,) = data
    N = length(ks[1])  # number of dimensions (usually 3)
    variance = zero(vn_rms)
    for i ∈ eachindex(ks, cs)
        k⃗ = SVector(ks[i])
        k² = sum(abs2, k⃗)
        @test k² > 0  # make sure zero mode is not included (even if we put kmin = 0.0)
        @test k⃗[1] ≥ 0  # Hermitia symmetry: modes kx < 0 are not included
        u⃗ = cs[i]
        u² = sum(abs2, u⃗)
        factor = iszero(k⃗[1]) ? 1 : 2  # Hermitian symmetry
        ku = norm(k⃗ × u⃗)               # just as a reference for comparison with k⃗ ⋅ u⃗ (which should be 0)
        @test k⃗ ⋅ u⃗ + ku ≈ ku          # check divergence-free (k⃗ ⋅ u⃗ = 0)
        variance += factor * u²
    end
    @test N * vn_rms^2 ≈ variance  # check variance
    nothing
end

@testset "Normal fluid forcing" begin
    rng = Xoshiro(42)
    forcing = @inferred ConstantFourierNormalFluidForcing(rng; α = 0.2f0, kmin = 0.0, kmax = 2.5)
    check_forcing(forcing)
end

