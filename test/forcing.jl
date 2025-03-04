using VortexPasta
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.PredefinedCurves
using VortexPasta.Forcing
using VortexPasta.SyntheticFields
using LinearAlgebra: ×
using Random
using Test

function generate_biot_savart_parameters(::Type{T}) where {T}
    Γ = 1.0
    a = 1e-6
    Δ = 1/4
    L = 2π
    Ls = (L, L, L)
    β = 3.0
    rcut = L / 2
    α = β / rcut
    kmax = 2α * β
    M = ceil(Int, kmax * L / π)
    Ns = (1, 1, 1) .* M
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short = NaiveShortRangeBackend(),
        backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4)),
        quadrature = GaussLegendre(3),
    )
end

function generate_filament(N)
    S = define_curve(TrefoilKnot(); translate = π, scale = π / 4)
    Filaments.init(S, ClosedFilament, N, QuinticSplineMethod())
end

# Check that the Fourier band forcing computes the right thing.
function check_fourier_band_forcing(forcing::FourierBandForcing, f::AbstractFilament, vL_computed, vs_self, cache_bs::BiotSavart.BiotSavartCache)
    vn_band = (forcing.vn)::FourierBandVectorField
    vs_band = (@inferred similar(vn_band))::FourierBandVectorField
    @test length(vs_band) == length(vn_band)  # same number of Fourier coefficients
    @test vs_band.qs == vn_band.qs    # same wavevectors
    @test vs_band.Δks == vn_band.Δks  # same wavevectors

    # Copy coefficients from vs_hat
    data_bs = @inferred BiotSavart.get_longrange_field_fourier(cache_bs)
    @test data_bs.state.quantity == :velocity
    @test data_bs.state.smoothed == true  # this means that the field needs to be unsmoothed (by reverting Gaussian filter)
    vs_hat = data_bs.field
    SyntheticFields.from_fourier_grid!(vs_band, vs_hat)  # copy coefficients of smoothed velocity field

    # Revert Gaussian filter
    α_ewald = cache_bs.params.α
    let (; cs, qs, Δks) = vs_band
        for i in eachindex(cs, qs)
            k⃗ = Δks .* qs[i]
            k² = sum(abs2, k⃗)
            cs[i] = cs[i] * exp(k² / (4 * α_ewald^2))  # deconvolve
        end
    end

    # Compute v_{ns} in Fourier band
    vns_band = similar(vn_band)
    @. vns_band.cs = vn_band.cs - vs_band.cs

    # Estimate vL
    vL = similar(vL_computed)
    for i in eachindex(vL, f)
        s⃗ = f[i]
        s⃗′ = f[i, UnitTangent()]
        v⃗ₙₛ = vns_band(s⃗)
        vL[i] = vs_self[i] + forcing.α * s⃗′ × v⃗ₙₛ
    end

    @test vL ≈ vL_computed

    nothing
end

@testset "Forcing" begin
    # Compute self-induced velocity of vortex filament
    N = 100
    f = generate_filament(N)
    fs = [f]
    p = generate_biot_savart_parameters(Float64)
    cache_bs = BiotSavart.init_cache(p, fs)
    vs = map(similar ∘ nodes, fs)
    velocity_on_nodes!(vs, cache_bs, fs)
    vs_self = map(copy, vs)
    α_ewald = p.α

    # Obtain coarse-grained velocity field in Fourier space
    data = BiotSavart.get_longrange_field_fourier(cache_bs)
    vs_hat = data.field
    state = data.state
    @test state.quantity == :velocity
    @test state.smoothed == true

    # Initialise random normal fluid velocity in Fourier space
    rng = Xoshiro(42)
    Ls = (4π, 2π, 2π)
    vn_rms = 1.0
    vn = @inferred FourierBandVectorField(undef, Ls; kmin = 0.1, kmax = 1.5)
    SyntheticFields.init_coefficients!(rng, vn, vn_rms)  # randomly set non-zero Fourier coefficients of the velocity field

    @testset "NormalFluidForcing" begin
        forcing = @inferred NormalFluidForcing(vn; α = 0.8, α′ = 0)
        cache = @inferred Forcing.init_cache(forcing, vs_hat)
        Forcing.update_cache!(cache, forcing, vs_hat, α_ewald)  # doesn't do anything for NormalFluidForcing
        @test startswith("NormalFluidForcing")(repr(forcing))
        for i in eachindex(fs, vs)
            copyto!(vs[i], vs_self[i])  # self-induced velocity (before forcing)
            Forcing.apply!(forcing, cache, vs[i], fs[i])
        end
    end

    @testset "FourierBandForcing" begin
        forcing = @inferred FourierBandForcing(vn; α = 0.8)
        cache = @inferred Forcing.init_cache(forcing, vs_hat)
        Forcing.update_cache!(cache, forcing, vs_hat, α_ewald)
        @test startswith("FourierBandForcing")(repr(forcing))
        for i in eachindex(fs, vs)
            copyto!(vs[i], vs_self[i])  # self-induced velocity (before forcing)
            Forcing.apply!(forcing, cache, vs[i], fs[i])
            check_fourier_band_forcing(forcing, fs[i], vs[i], vs_self[i], cache_bs)
        end
    end
end
