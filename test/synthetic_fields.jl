using VortexPasta
using VortexPasta.SyntheticFields
using StaticArrays
using StructArrays
using LinearAlgebra
using FFTW
using Random
using Test

# Note: this may fail or give wrong results if ûs is not large enough.
function to_fourier_grid!(ûs::StructArray, field::FourierBandVectorField{T, N}) where {T, N}
    SyntheticFields.to_fourier_grid!(StructArrays.components(ûs), field)
    ûs
end

function check_synthetic_field(field::FourierBandVectorField{T, N}; u_rms, kmin, kmax) where {T, N}
    (; cs, qs, Δks,) = field
    variance = zero(T)
    divergence = zero(T)
    k²_min = T(Inf)
    kx_min = T(Inf)
    k²_max = zero(T)
    for i ∈ eachindex(qs, cs)
        q⃗ = SVector(qs[i])
        k⃗ = (q⃗ .* Δks)::SVector{N, T}
        k² = sum(abs2, k⃗)
        k²_min = min(k²_min, k²)
        k²_max = max(k²_max, k²)
        kx_min = min(kx_min, k⃗[1])
        u⃗ = cs[i]
        u² = sum(abs2, u⃗)
        factor = 1 + SyntheticFields.has_implicit_conjugate(q⃗)  # Hermitian symmetry
        divergence += abs2(k⃗ ⋅ u⃗)
        variance += factor * u²
    end
    @test kx_min ≥ 0  # Hermitian symmetry: modes kx < 0 are not included
    @test k²_min ≥ kmin^2
    @test k²_max ≤ kmax^2
    @test divergence < eps(variance)^2  # divergence is zero
    @test N * u_rms^2 ≈ variance        # check variance
    nothing
end

function verify_evaluation(field::FourierBandVectorField{T, N}, Ls::NTuple{N, T}; u_rms, kmax) where {T, N}
    # Directly evaluate field on a grid
    Ns = map(Ls) do L
        Δk_half = T(π) / L
        floor(Int, kmax / Δk_half + 1)  # enough Fourier modes to include the forcing band
    end
    Δxs = Ls ./ Ns

    let S = T === Float64 ? Float32 : Float64  # select a different type as input
        x⃗ = zero(SVector{N, S})
        @test @inferred(field(x⃗)) isa SVector{N, T}  # check that the returned type is T and not S
    end

    us = StructArray{SVector{N, T}}(undef, Ns)
    u⃗_mean = zero(eltype(us))
    u⃗_var = zero(eltype(us))
    for I ∈ CartesianIndices(us)
        x⃗ = Tuple(I - oneunit(I)) .* Δxs
        u⃗ = field(x⃗)
        u⃗_mean = @. u⃗_mean + u⃗
        u⃗_var = @. u⃗_var + u⃗ * u⃗
        us[I] = u⃗
    end
    Nprod = length(us)
    u⃗_mean = @. u⃗_mean / Nprod
    u⃗_var = @. u⃗_var / Nprod
    u⃗_rms = @. sqrt(u⃗_var)
    @test u⃗_mean + u⃗_rms ≈ u⃗_rms rtol=eps(T)  # the mean is practically zero (assumes kmin > 0)
    @test sum(u⃗_var) / N ≈ u_rms^2

    # Perform FFT and compare with field obtained using SyntheticFields.to_fourier_grid!.
    ûs_from_grid = StructArray{SVector{N, Complex{T}}}(rfft.(StructArrays.components(us))) ./ Nprod
    Ms = ntuple(d -> d == 1 ? Ns[d] ÷ 2 + 1 : Ns[d], Val(N))
    ûs = StructArray{SVector{N, Complex{T}}}(undef, Ms)
    to_fourier_grid!(ûs, field)

    @test ûs ≈ ûs_from_grid

    # Check variance in Fourier space
    u⃗_var_fft = zero(u⃗_var)
    for I ∈ CartesianIndices(ûs)
        factor = isone(I[1]) ? 1 : 2
        û = ûs[I]
        u⃗_var_fft = @. u⃗_var_fft + factor * abs2(û)
    end
    @test sum(u⃗_var_fft) / N ≈ u_rms^2

    nothing
end

function test_io(field::FourierBandVectorField)
    mktempdir() do path
        (; cs, qs,) = field
        cs_prev = copy(cs)
        qs_prev = copy(qs)
        fname = joinpath(path, "field.h5")
        SyntheticFields.save_coefficients(fname, field)
        empty!(cs)  # they will be resized in load_coefficients!
        empty!(qs)
        SyntheticFields.load_coefficients!(fname, field)
        @test cs == cs_prev
        @test qs == qs_prev
    end
end

@testset "Synthetic fields ($T)" for T ∈ (Float32, Float64)
    @testset "FourierBandVectorField" begin
        rng = Xoshiro(42)
        Ls = T.((2π, 2π, 8π))  # elongated domain
        u_rms = 2.0; kmin = 0.1; kmax = 1.5;
        field = @inferred FourierBandVectorField(undef, Ls; kmin, kmax)
        N = length(field)
        SyntheticFields.init_coefficients!(rng, field, u_rms)
        @test N == length(field.cs) == length(field.qs)
        @test startswith("FourierBandVectorField{$T, 3}")(repr(field))
        check_synthetic_field(field; u_rms, kmin, kmax)
        verify_evaluation(field, Ls; u_rms, kmax)
        test_io(field)
        @testset "Remove coefficients" begin
            nremoved = SyntheticFields.remove_zeros!(field)
            @test nremoved == 0  # no coefficients are currently zero
            @test length(field) == N  # no coefficients have been removed
            # Now set a coefficient to zero
            field.cs[end - 2] = zero(eltype(field.cs))
            nremoved = SyntheticFields.remove_zeros!(field)
            @test nremoved == 1
            @test length(field) == N - 1
        end
    end
    @testset "UniformVectorField" begin
        u⃗ = T.((2.0, 1.2, -2.4))
        field = @inferred UniformVectorField(u⃗)
        @test startswith("UniformVectorField{$T, 3} with value")(repr(field))
        x⃗ = T.((0.2, 0.5, 1.5))
        @test (@inferred field(x⃗)) === SVector(u⃗)
    end
end
