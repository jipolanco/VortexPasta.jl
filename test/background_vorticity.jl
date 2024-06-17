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
using VortexPasta.BasicTypes: VectorOfVectors
using LinearAlgebra: norm
using Test

function evaluate_bs_on_nodes!(fields, fs; L, α, β)
    rcut = β / α
    kmax = 2α * β
    M = floor(Int, L * kmax / π)
    params = ParamsBiotSavart(
        Float64;
        Γ = 1.0, a = 1e-8,
        α, Ls = L,
        rcut, Ns = (M, M, M),
        backend_long = NonuniformFFTsBackend(m = HalfSupport(4))  # ~1e-6 accuracy
    )
    cache = BiotSavart.init_cache(params, fs)
    BiotSavart.compute_on_nodes!(fields, cache, fs)
    fields
end

##

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

αs = [32, 16] ./ L

fields_comp = map(αs) do α
    β = 3.5  # ~ 1e-6 accuracy
    evaluate_bs_on_nodes!(fields, fs; L, α, β)
    deepcopy(fields)
end

# @show norm(fields_comp[1].velocity - fields_comp[2].velocity) / norm(fields_comp[1].velocity)
# @show norm(fields_comp[1].streamfunction - fields_comp[2].streamfunction) / norm(fields_comp[1].streamfunction)

@testset "Independence on α" begin
    @test isapprox(fields_comp[1].velocity, fields_comp[2].velocity; rtol = 5e-7)
    @test isapprox(fields_comp[1].streamfunction, fields_comp[2].streamfunction; rtol = 2e-8)
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
