using Test
using VortexPasta.BiotSavart
using .BiotSavart: ChebyshevApproximations
using SpecialFunctions: erf, erfc

# Approximate a Gaussian up to some cut-off distance rcut.
function test_gaussian(::Type{T}; symmetry = :none) where {T}
    β = T(6.0)
    rcut = T(0.4)
    α = β / rcut
    rtol = 10 * eps(T)
    f = ChebyshevApproximations.approximate(rcut; symmetry = Val(symmetry), rtol) do r
        2 * α / sqrt(T(π)) * exp(-α^2 * r^2)  # the factor 2 is such that ∫_0^∞ f(r) dr = 1
    end
    F = ChebyshevApproximations.integrate(f)

    # Verify that the integral corresponds to erf.
    rs = range(0, rcut; length = 1001)
    @test F.(rs) ≈ erf.(α .* rs) rtol=rtol

    # Verify that ∫_0^∞ (r²/2) * f(r) dr = 1 / (4 * α²) [used in background vorticity correction]
    # Note that we truncate the integral at rcut, which leads to a small truncation error (very small when β = 6).
    # 1. Using a new Chebyshev approximation
    f_rr = ChebyshevApproximations.approximate(r -> r^2 * f(r) / 2, rcut; symmetry = Val(symmetry), rtol)
    F_rr = ChebyshevApproximations.integrate(f_rr)
    F_rr_integral_expected = 1 / (4 * α^2)
    @test F_rr(rcut) ≈ F_rr_integral_expected rtol=2rtol

    # 2. Using multiplication of Chebyshev polynomials to integrate x * F(x).
    # Should be faster but seems to be slightly less accurate (10x).
    let cs = F.cs
        if symmetry == :none
            @assert F isa ChebyshevApproximations.ChebyshevSeries{:none}
            v = cs[2] / 2  # = c₁/4 * ∫_{-1}^{1} T₀(x) dx
            for i in eachindex(cs)[4:2:end]  # cs[i] = c_{i - 1}
                n = i - 2
                v += (cs[i - 2] + cs[i]) / (2 * (1 - n^2))
            end
        elseif symmetry == :even
            # In this case F has odd symmetry, and cs only contains coefficients c_{2i - 1}
            @assert F isa ChebyshevApproximations.ChebyshevSeries{:odd}
            v = cs[1] / 2
            for i in eachindex(cs)[2:end]  # cs[i] = c_{2i - 1}
                n = 2i - 2
                v += (cs[i - 1] + cs[i]) / (2 * (1 - n^2))
            end
        end
        integ = rcut^2 * (T(1/2) - v)
        @test integ ≈ F_rr_integral_expected rtol=20rtol
    end

    nothing
end

@testset "ChebyshevApproximations" begin
    @testset "Precision: $T" for T in (Float32, Float64)
        @testset "Symmetry: $symmetry" for symmetry in (:none, :even)
            test_gaussian(T; symmetry)
        end
    end
end
