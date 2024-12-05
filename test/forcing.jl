using VortexPasta.Forcing
using VortexPasta.SyntheticFields
using Random
using Test

@testset "NormalFluidForcing" begin
    rng = Xoshiro(42)
    Ls = (4π, 2π, 2π)
    vn_rms = 1.0
    vn = FourierBandVectorField(undef, Ls; kmin = 0.1, kmax = 1.5)
    SyntheticFields.init_coefficients!(rng, vn, vn_rms)  # randomly set non-zero Fourier coefficients of the velocity field
    forcing = NormalFluidForcing(vn; α = 0.8, α′ = 0)
    @test startswith("NormalFluidForcing")(repr(forcing))
end
