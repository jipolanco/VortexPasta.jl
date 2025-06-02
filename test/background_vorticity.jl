# Test case with non-zero mean vorticity associated to vortex filaments.
# This non-zero mean vorticity is compensated by a background vorticity.
#
# In the long-range part of Ewald summation, the background vorticity is implicitly included
# by setting ω̂(k = 0) = 0.
#
# In the short-range part, there is nothing to be done when computing the velocity.
# Extra care is needed when computing the streamfunction or the kinetic energy, so that
# results don't depend on the splitting parameter α.

using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.PredefinedCurves
using VortexPasta.Containers: VectorOfVectors
using LinearAlgebra: norm
using Test

function init_params_biot_savart(; L, Ngrid, β)
    kmax = Ngrid * π / L
    α = kmax / (2β)
    rcut = β / α
    ParamsBiotSavart(
        Float64;
        Γ = 1.0, a = 1e-8,
        α, Ls = L,
        rcut, Ns = (Ngrid, Ngrid, Ngrid),
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4))  # ~1e-6 accuracy
    )
end

function evaluate_bs_on_nodes!(fields, fs, params)
    cache = BiotSavart.init_cache(params, fs)
    BiotSavart.compute_on_nodes!(fields, cache, fs)
    cache
end

## Check independence of results on Ewald's splitting parameter α

Np = 64
L = 2π

# Create 4×4 array of aligned vortices with identical helical perturbations.
Nx = 4
xs = range(0, L; length = 2 * Nx + 1)[2:2:end]
ys = xs

fs = map(Iterators.product(xs, ys)) do (x, y)
    p = PeriodicLine(r = t -> cispi(4t) / 100)
    S = define_curve(p; scale = L, translate = (x, y, L/2), orientation = +1)
    ζs = range(0, 1; length = 2 * Np + 1)[2:2:end]
    Filaments.init(S, ClosedFilament, ζs, QuinticSplineMethod())
end |> vec

vs = map(similar ∘ nodes, fs) |> VectorOfVectors
ψs = map(similar ∘ nodes, fs) |> VectorOfVectors
fields = (velocity = vs, streamfunction = ψs)

Ns = [32, 64]  # this also changes the splitting parameter α by a factor 2

fields_comp = map(Ns) do Ngrid
    β = 3.5  # ~ 1e-6 accuracy
    params = init_params_biot_savart(; L, Ngrid, β)
    evaluate_bs_on_nodes!(fields, fs, params)
    deepcopy(fields)
end

# @show norm(fields_comp[1].velocity - fields_comp[2].velocity) / norm(fields_comp[1].velocity)
# @show norm(fields_comp[1].streamfunction - fields_comp[2].streamfunction) / norm(fields_comp[1].streamfunction)

@testset "Independence on α" begin
    @test isapprox(fields_comp[1].velocity, fields_comp[2].velocity; rtol = 5e-7)
    @test isapprox(fields_comp[1].streamfunction, fields_comp[2].streamfunction; rtol = 5e-7)
end

# using GLMakie
#
# fig = Figure()
# ax = Axis3(fig[1, 1]; aspect = :data)
# hidespines!(ax)
# wireframe!(ax, Rect(0, 0, 0, L, L, L); color = (:grey, 0.5))
# for f ∈ fs
#     plot!(ax, f; markersize = 4)
# end
# fig

## Single vortex + Gaussian smoothing

S = define_curve(PeriodicLine(); scale = L, translate = (L/2, L/2, L/2))
f = Filaments.init(S, ClosedFilament, Np, CubicSplineMethod())
fs = [f]

vs = map(similar ∘ nodes, fs) |> VectorOfVectors
ψs = map(similar ∘ nodes, fs) |> VectorOfVectors
fields = (velocity = vs, streamfunction = ψs)
β = 3.5  # ~ 1e-6 accuracy
Ngrid = 32
params = init_params_biot_savart(; L, Ngrid, β)
cache = evaluate_bs_on_nodes!(fields, fs, params)

@testset "Coarse-grained vorticity" begin
    # Now that we have computed the actual fields in Fourier space, obtain the coarse-grained
    # vorticity at a chosen scale.
    ℓ = L / 8
    (; field, wavenumbers, state,) = BiotSavart.to_coarse_grained_vorticity!(cache.longrange, ℓ)

    @test state.quantity == :vorticity
    @test state.smoothing_scale == ℓ

    ω_back = -params.Γ / L^2  # background vorticity required to neutralise the vortex charge

    # Expected value of coarse-grained vorticity evaluated on the vortex (including
    # background vorticity). The first term simply corresponds to a 2D Gaussian kernel
    # (evaluated at the centre), since we're dealing with a straight vortex.
    ω_expected = params.Γ / (2π * ℓ^2) + ω_back

    # Now interpolate values on filament
    ωs_ℓ = similar(vs)
    BiotSavart.interpolate_to_physical!(cache.longrange)     # interpolate to vortex positions
    BiotSavart.copy_long_range_output!(ωs_ℓ, cache.longrange)  # copy results

    @test all(ω -> isapprox(ω, ωs_ℓ[1][1]; rtol = 1e-6), ωs_ℓ[1])  # all vorticities are equal (up to chosen accuracy)
    ω⃗ = sum(ωs_ℓ[1]) ./ length(ωs_ℓ[1])  # average vorticity
    @test ω⃗[3] ≈ ω_expected rtol=1e-5
end
