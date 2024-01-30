using Test
using VortexPasta.Filaments
using VortexPasta.PredefinedCurves
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using Rotations: Rotations

function generate_biot_savart_parameters(::Type{T}) where {T}
    # Parameters are relevant for HeII if we interpret dimensions in cm and s.
    Γ = 9.97e-4
    a = 1e-8
    Δ = 1/4
    L = 2π
    Ngrid = 64
    Ls = (L, L, L)
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax = (Ngrid ÷ 2) * 2π / L
    α::T = kmax / 8
    rcut = 4 / α
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short = CellListsBackend(2),
        backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4)),
        quadrature = GaussLegendre(2),
    )
end

function generate_linked_rings(
        ::Type{T}, N, method = QuinticSplineMethod();
        R,
    ) where {T}
    p = Ring()
    δ = R/2
    scale = R
    rotate = Rotations.RotX(π/2)
    Sl = define_curve(p; scale, translate = (-δ, 0, 0))           # left ring
    Sr = define_curve(p; scale, rotate, translate = (+δ, 0, 0))  # right ring
    [
        Filaments.init(Sl, ClosedFilament{T}, N, method),
        Filaments.init(Sr, ClosedFilament{T}, N, method),
    ]
end

function test_linked_rings(::Type{T}, N, method; R) where {T}
    fs = @inferred generate_linked_rings(T, N, method; R)
    params = @inferred generate_biot_savart_parameters(T)

    tspan = (0.0, 3.2)  # not used
    prob = VortexFilamentProblem(fs, tspan, params)
    iter = init(prob, RK4(); dt = 0.1)

    (; Γ,) = params
    E = @inferred Diagnostics.kinetic_energy(iter; quad = GaussLegendre(2))
    H_no_quad = @inferred Diagnostics.helicity(iter; quad = nothing)
    H = @inferred Diagnostics.helicity(iter; quad = GaussLegendre(2))
    @test E isa T
    @test H_no_quad isa T
    @test H isa T

    # @show (H - H_no_quad) / H
    @test isapprox(H, H_no_quad)  # there is really no difference between the two

    # In the case of two unknotted linked rings, the helicity is H = 2 Γ² L where L is the
    # linking number. For two linked rings, |L| = 1.
    linking = H / (2 * Γ^2)

    # Note: the accuracy seems to be mainly controlled by the splitting parameter α/kmax
    # and by the accuracy of the NUFFTs.
    # @show linking + 1
    rtol = 2e-6
    @test isapprox(linking, -1; rtol)

    nothing
end

@testset "Linked rings" begin
    T = Float32
    N = 64
    method = QuinticSplineMethod()
    R = 1.2
    test_linked_rings(T, N, method; R)
end
