using Test
using Random
using StableRNGs: StableRNG
using VortexPasta
using VortexPasta.BiotSavart
using LinearAlgebra: norm
using Adapt: adapt
using KernelAbstractions: KernelAbstractions as KA, POCLBackend

function generate_point_vortices(rng, N, Ls::NTuple{2, T}; neutral = false) where {T <: AbstractFloat}
    if neutral && isodd(N)
        error("number of vortices must be even for neutral charge")
    end
    points = Vector{Vec2{T}}(undef, N)
    signs = similar(points, SignValue)
    for j in eachindex(points)
        points[j] = rand(rng, Vec2{T}) .* Ls  # rescale positions
        signs[j] = rand(rng, Bool) ? UP : DOWN  # randomly choose UP or DOWN
    end
    if neutral
        signs[begin:(end ÷ 2)] .= UP
        signs[(end ÷ 2 + 1):end] .= DOWN
    end
    PointVortexState(points, signs)
end

function point_vortex_interactions(vortices::PointVortexState, params::ParamsPointVortex{T}; kws...) where {T}
    cache = @inferred BiotSavart.init_cache(params)
    fields = (
        velocity = similar(vortices.points),
        streamfunction = similar(vortices.points, T),
    )
    compute_on_nodes!(fields, cache, vortices; kws...)
    energy = kinetic_energy_2D(fields.streamfunction, vortices.signs, params)
    (; cache, params, vortices, energy, fields...)
end

# Kinetic energy per unit mass
function kinetic_energy_2D(ψs, signs, params)
    (; Ls, Γ) = params
    T = typeof(Γ)
    E = zero(T)
    for i in eachindex(ψs, signs)
        E += T(BiotSavart.to_real(signs[i]) * ψs[i])
    end
    E = E * T(Γ / 2 / prod(Ls))
    E
end

##

T = Float64
Lx, Ly = T(2π), T(3π)
Ls = (Lx, Ly)
Γ = 2.0

N = 21  # odd number to make sure the charge is not neutral
rng = StableRNG(42)
vortices = generate_point_vortices(rng, N, Ls)

# Base run
backend = CPU()
backend_short = CellListsBackend(backend)
backend_long = NonuniformFFTsBackend(backend; σ = 1.5, m = HalfSupport(4))
splitting = @inferred GaussianSplitting(; Ls, rcut = Lx / 4, β = 4.0)
params = @inferred ParamsPointVortex(; Γ, splitting, backend_short, backend_long)
ref = point_vortex_interactions(vortices, params)

@testset "Computing only the velocity" begin
    vs_only = @inferred velocity_on_nodes(ref.cache, ref.vortices)
    @test vs_only ≈ ref.velocity rtol=1e-15
end

@testset "callback_vorticity" begin
    E_from_vorticity = Ref(zero(T))
    function callback_vorticity(cache::LongRangeCache{2})
        (; field, wavenumbers, state) = BiotSavart.get_longrange_field_fourier(cache)
        @assert state.quantity == :vorticity
        @assert state.smoothing_scale == 0  # unsmoothed field
        # To make things simple, we copy data to the CPU if it's on the GPU.
        wavenumbers = adapt(Array, wavenumbers)
        # Note: in 2D, vorticity is contained in field[1] (field[2] is unused)
        ω_hat = adapt(Array, field[1])::AbstractArray{<:Complex, 2}  # ωz(kx, ky) in Fourier space
        with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)  # this depends on the long-range backend
        @assert with_hermitian_symmetry == (wavenumbers[1][end] > 0)
        E = 0.0
        for I in CartesianIndices(ω_hat)
            k⃗ = map(getindex, wavenumbers, Tuple(I))
            kx = k⃗[1]
            factor = (!with_hermitian_symmetry || kx == 0) ? 0.5 : 1.0
            k² = sum(abs2, k⃗)
            if !iszero(k²)
                E += factor * abs2(ω_hat[I]) / k²
            end
        end
        E_from_vorticity[] = E  # update value of "global" variable
        nothing
    end
    point_vortex_interactions(vortices, params; callback_vorticity)
    @test E_from_vorticity[] > 0
end

@testset "Float32" begin
    points_f32 = Vec2{Float32}.(ref.vortices.points)
    vortices_f32 = PointVortexState(points_f32, ref.vortices.signs)
    @test vortices_f32 isa PointVortexState{Float32}
    params_f32 = @inferred ParamsPointVortex(; Γ = Float32(Γ), splitting, backend_short, backend_long)
    @test params isa ParamsPointVortex{Float32}
    results = point_vortex_interactions(vortices_f32, params_f32)
    vs_only = @inferred velocity_on_nodes(results.cache, vortices_f32)
    @test typeof(vs_only) == typeof(vortices_f32.points)  # returned velocities are Float32
    # @show norm(results.velocity - ref.velocity) / norm(ref.velocity)
    # @show norm(results.streamfunction - ref.streamfunction) / norm(ref.streamfunction)
    # @show (results.energy - ref.energy) / ref.energy
    @test results.energy ≈ ref.energy rtol=1f-5
    @test results.velocity ≈ ref.velocity rtol=1f-5
    @test results.streamfunction ≈ ref.streamfunction rtol=1f-5
end

@testset "POCLBackend" begin
    backend = POCLBackend()
    backend_short = CellListsBackend(backend, 2)
    backend_long = NonuniformFFTsBackend(backend; σ = 1.5, m = HalfSupport(4), gpu_method = :global_memory)
    params = @inferred ParamsPointVortex(; Γ, splitting, backend_short, backend_long)
    results = point_vortex_interactions(ref.vortices, params)
    @test KA.get_backend(results.cache.longrange) isa POCLBackend
    @test KA.get_backend(results.cache.shortrange) isa POCLBackend
    # @show (results.energy - ref.energy) / ref.energy
    # @show norm(results.velocity - ref.velocity) / norm(ref.velocity)
    # @show norm(results.streamfunction - ref.streamfunction) / norm(ref.streamfunction)
    @test results.energy ≈ ref.energy rtol=1e-15
    @test results.velocity ≈ ref.velocity rtol=1e-15
    @test results.streamfunction ≈ ref.streamfunction rtol=1e-15
end

@testset "backend_long = ExactSumBackend()" begin
    backend_long = ExactSumBackend()
    params = @inferred ParamsPointVortex(; Γ, splitting, backend_short, backend_long)
    results = point_vortex_interactions(ref.vortices, params)
    # @show (results.energy - ref.energy) / ref.energy
    # @show norm(results.velocity - ref.velocity) / norm(ref.velocity)
    # @show norm(results.streamfunction - ref.streamfunction) / norm(ref.streamfunction)
    @test results.energy ≈ ref.energy rtol=1e-6
    @test results.velocity ≈ ref.velocity rtol=1e-7
    @test results.streamfunction ≈ ref.streamfunction rtol=1e-6
end

@testset "backend_short = NaiveShortRangeBackend()" begin
    backend_short = NaiveShortRangeBackend()
    params = @inferred ParamsPointVortex(; Γ, splitting, backend_short, backend_long)
    results = point_vortex_interactions(ref.vortices, params)
    # @show (results.energy - ref.energy) / ref.energy
    # @show norm(results.velocity - ref.velocity) / norm(ref.velocity)
    # @show norm(results.streamfunction - ref.streamfunction) / norm(ref.streamfunction)
    @test results.energy ≈ ref.energy rtol=1e-16
    @test results.velocity ≈ ref.velocity rtol=1e-16
    @test results.streamfunction ≈ ref.streamfunction rtol=1e-16
end

@testset "Results independent of rcut" begin
    rcuts = (Lx / 4, Lx / 6)
    A, B = map(rcuts) do rcut
        splitting = GaussianSplitting(; Ls, rcut, β = 4.0)
        params = ParamsPointVortex(; Γ = 2.0, splitting, backend_short, backend_long)
        point_vortex_interactions(ref.vortices, params)
    end
    # @show norm(B.velocity - A.velocity) / norm(A.velocity)
    # @show norm(B.streamfunction - A.streamfunction) / norm(A.streamfunction)
    # @show (A.energy - B.energy) / B.energy
    @test A.energy ≈ B.energy rtol=1e-7
    @test A.velocity ≈ B.velocity rtol=1e-7
    @test A.streamfunction ≈ B.streamfunction rtol=1e-6
end
