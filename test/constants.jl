# Test Zero and Infinity types.

using VortexPasta.Constants: Zero, Infinity
using Test

@testset "Zero and Infinity" begin
    @testset "Conversion" begin
        # Convert to same type
        @test convert(Zero, Zero()) === Zero()
        @test convert(Infinity, Infinity()) === Infinity()

        # Convert to float types
        for T ∈ (Float32, Float64)
            @test convert(T, Zero()) === zero(T)
            @test convert(T, Infinity()) === T(Inf)
        end

        # Convert Zero to integer types
        for T ∈ (Int32, Int64, UInt32)
            @test convert(T, Zero()) === zero(T)
        end
    end

    @testset "Product" begin
        for x ∈ (Zero(), 0.3, 100)
            @test Zero() === @inferred x * Zero()
            @test Zero() === @inferred Zero() * x
        end
        for x ∈ (Infinity(), 0.3, 100)
            @test Infinity() === @inferred x * Infinity()
            @test Infinity() === @inferred Infinity() * x
        end
    end

    @testset "Division" begin
        for x ∈ (0.3, 100)
            @test @inferred(x / Zero()) === Infinity()
            @test @inferred(Zero() / x) === Zero()
            @test @inferred(Zero() ÷ x) === Zero()
            @test @inferred(x / Infinity()) === Zero()
            @test @inferred(Infinity() / x) === Infinity()
            @test @inferred(Infinity() ÷ x) === Infinity()
        end
    end

    @testset "Sum / subtraction" begin
        for x ∈ (Infinity(), 0.3, 100)
            @test Infinity() === @inferred x + Infinity()
            @test Infinity() === @inferred Infinity() + x
            @test Infinity() === @inferred Infinity() - x
        end
    end

    @testset "Comparison (<, >, ...)" begin
        @test (Infinity() < Infinity()) === false
        @test (Infinity() > Infinity()) === false
        for x ∈ (0.3, 100)
            @test true === @inferred Infinity() > x
            @test true === @inferred x < Infinity()
        end
    end

    @testset "Functions" begin
        @test 1 == @inferred exp(Zero())
    end
end
