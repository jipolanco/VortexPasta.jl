# Verify correctness of SIMDFunctions implementations
using VortexPasta.BiotSavart: SIMDFunctions
using SpecialFunctions: SpecialFunctions
using Test

function test_exp(xs)
    T = eltype(xs)
    rtol = 1 * eps(T)  # very accurate! (1 ulp)
    for x in xs
        @test exp(x) ≈ SIMDFunctions.exp(x) rtol=rtol
    end
    nothing
end

function test_erf(xs)
    T = eltype(xs)
    rtol = 4 * eps(T)
    for x in xs
        @test SpecialFunctions.erf(x) ≈ SIMDFunctions.erf(x) rtol=rtol
    end
    nothing
end

function test_erfc(xs)
    T = eltype(xs)
    # Note: erfc computed as 1 - erf gets quite inaccurate (in a relative sense) when we
    # reach values very close to 0. But that's ok, since we only care about absolute error.
    atol = 4 * eps(T)
    for x in xs
        @test SpecialFunctions.erfc(x) ≈ (one(x) - SIMDFunctions.erf(x)) atol=atol
    end
    nothing
end

@testset "SIMDFunctions" begin
    @testset "exp" begin
        @testset "Float64" begin
            test_exp(range(-200.0, 200.0; length = 1001))
            test_exp(logrange(1e-12, 1e12; length = 1001))
            test_exp(-logrange(1e-12, 1e12; length = 1001))
        end
        @testset "Float32" begin
            test_exp(range(-200f0, 200f0; length = 1001))
            test_exp(logrange(1f-12, 1f12; length = 1001))
            test_exp(-logrange(1f-12, 1f12; length = 1001))
        end
    end
    @testset "erf" begin
        @testset "Float64" test_erf(range(0.0, 6.0; length = 1001))      # erf(6.0) == 1.0
        @testset "Float32" test_erf(range(0.0f0, 4.0f0; length = 1001))  # erf(4.0f0) == 1.0f0
    end
    @testset "erfc" begin
        @testset "Float64" test_erfc(range(0.0, 6.0; length = 1001))      # erf(6.0) == 1.0
        @testset "Float32" test_erfc(range(0.0f0, 4.0f0; length = 1001))  # erf(4.0f0) == 1.0f0
    end
end
