using Test
using Random
using StableRNGs: StableRNG
using VortexPasta
using VortexPasta.BiotSavart
using LinearAlgebra: norm
using KernelAbstractions: POCLBackend

# TODO: test:
# - different short-range and long-range backends
# - independence on Ewald parameters
# - POCLBackend
# - callback_vorticity (e.g. to get energy spectrum)
# - computing only velocity

N = 21  # odd number to make sure the charge is not neutral
T = Float64
Lx, Ly = T(2π), T(3π)
Ls = (Lx, Ly)
rng = StableRNG(42)
vortices = Vector{Vec2{T}}(undef, N)
signs = similar(vortices, SignValue)
for j in eachindex(vortices)
    vortices[j] = rand(rng, Vec2{T}) .* Ls  # rescale positions
    signs[j] = rand(rng, Bool) ? UP : DOWN  # randomly choose UP or DOWN
end

# Ensure neutral charge
# if isodd(length(vortices))
#     pop!(vortices)
#     pop!(signs)
#     signs[begin:(end ÷ 2)] .= UP
#     signs[(end ÷ 2 + 1):end] .= DOWN
#     N = length(vortices)
# end

npos = count(==(UP), signs)
nneg = length(signs) - npos

backend = CPU()
backend_short = NaiveShortRangeBackend(backend)
backend_long = NonuniformFFTsBackend(backend)
# backend_long = ExactSumBackend()
splitting = @inferred GaussianSplitting(; Ls, rcut = Lx / 4, β = 4.0)
params = @inferred ParamsPointVortex(; Γ = 2.0, splitting, backend_short, backend_long)

cache = @inferred BiotSavart.init_cache(params)

fields = (
    velocity = similar(vortices),
    streamfunction = similar(vortices, T),
)

compute_on_nodes!(fields, cache, vortices, signs)

@testset "Computing only the velocity" begin
    vs_only = @inferred velocity_on_nodes(cache, vortices, signs)
    @test vs_only ≈ fields.velocity rtol=1e-15
end

@testset "Results independent of rcut" begin
    rcuts = (Lx / 4, Lx / 6)
    A, B = map(rcuts) do rcut
        splitting = @inferred GaussianSplitting(; Ls, rcut, β = 4.0)
        params = @inferred ParamsPointVortex(; Γ = 2.0, splitting, backend_short, backend_long)
        cache = @inferred BiotSavart.init_cache(params)
        fields = (
            velocity = similar(vortices),
            streamfunction = similar(vortices, T),
        )
        compute_on_nodes!(fields, cache, vortices, signs)
        (; fields..., splitting)
    end
    # @show norm(B.velocity - A.velocity) / norm(A.velocity)
    # @show norm(B.streamfunction - A.streamfunction) / norm(A.streamfunction)
    @test A.velocity ≈ B.velocity rtol=1e-6
    @test A.streamfunction ≈ B.streamfunction rtol=1e-6
end

# pointdata = BiotSavart.PointVortexData(T)
# BiotSavart.init_cache_short(params, pointdata)
# longrange = @inferred BiotSavart.init_cache_long(params, pointdata)
