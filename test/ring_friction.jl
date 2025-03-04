# Test vortex ring affected by mutual friction due to a quiescent normal fluid.

using Test
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.PredefinedCurves
using VortexPasta.BiotSavart
using VortexPasta.SyntheticFields: UniformVectorField
using VortexPasta.Forcing: Forcing, NormalFluidForcing
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using LinearAlgebra: norm, normalize, ⋅, ×
using UnicodePlots: UnicodePlots, lineplot, lineplot!

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

function test_ring_friction_static(f, params, forcing::NormalFluidForcing)
    T = eltype(params)
    @test T <: AbstractFloat
    vs = similar(nodes(f))  # velocities
    @test eltype(f) === Vec3{T}
    @test eltype(vs) === Vec3{T}
    fs = [f]
    vs_all = [vs]
    cache = @inferred BiotSavart.init_cache(params, fs)

    # Compute self-induced velocity
    BiotSavart.velocity_on_nodes!(vs_all, cache, fs)

    # Check that velocities are all in the +Z direction (ring translation velocity).
    vs_orig = copy(vs)  # velocities before mutual friction
    vz_orig = getindex.(vs_orig, 3)
    @test norm.(vs_orig) ≈ vz_orig rtol=eps(T)

    # Now modify velocities due to quiescent normal fluid
    (; vn, α, α′,) = forcing
    forcing_cache = @inferred Forcing.init_cache(forcing, cache)
    @test iszero(vn.(nodes(f)))     # field is zero everywhere (in particular at vortex points)
    Forcing.apply!(forcing, forcing_cache, vs, f)  # modify vortex velocities

    # Subtract original velocity to obtain just the one due to mutual friction.
    vf = vs - vs_orig

    # Obtain radial component (due to Magnus force) and the rest (due to drag)
    vf_r = similar(vz_orig)  # radial component of vf
    vf_rest = similar(vf)    # = vf minus its radial component
    for i ∈ eachindex(vf, f)
        ρ⃗ = f[i, CurvatureVector()]
        r̂ = -normalize(ρ⃗)  # radial direction (points outwards)
        vr = vf[i] ⋅ r̂     # radial forcing velocity
        vf_r[i] = vr
        vf_rest[i] = vf[i] - vr * r̂
    end

    # The radial velocity (due to Magnus) points inwards, tending to shrink the vortex.
    # It's proportional to the α coefficient.
    @test vf_r ≈ -α * vz_orig

    # The drag force should be in the -Z direction (opposite to the ring translation).
    vf_z = getindex.(vf_rest, 3)
    @test norm.(vf_rest) ≈ -vf_z                             # vf_rest points in the -Z direction
    @test vf_z ≈ -α′ * vz_orig rtol=max(sqrt(eps(T)), 1e-7)  # it's proportional to the α′ coefficient

    nothing
end

function test_ring_friction_dynamic(f, params, forcing::NormalFluidForcing)
    T = eltype(params)
    tspan = (T(0), T(0.1))
    prob = VortexFilamentProblem([f], tspan, params)

    δ = minimum(node_distance, prob.fs)
    refinement = RefineBasedOnSegmentLength(0.75 * δ)
    dt = BiotSavart.kelvin_wave_period(params, δ) / 2
    scheme = RK4()

    time = T[]
    length = T[]
    energy = T[]
    impulse = Vec3{T}[]

    function callback(iter)
        local (; fs, vs, vn, vL, t, nstep,) = iter
        local quad = GaussLegendre(4)
        local L = Diagnostics.filament_length(iter; quad)
        local E = Diagnostics.kinetic_energy(iter; quad)
        local p⃗ = Diagnostics.vortex_impulse(iter; quad)
        push!(time, t)
        push!(length, L)
        push!(energy, E)
        push!(impulse, p⃗)
        nothing
    end

    iter = @inferred init(
        prob, scheme;
        dt, refinement,
        forcing,
        adaptivity = NoAdaptivity(),
        callback,
    )

    @test contains("NormalFluidForcing")(repr(iter))

    solve!(iter)

    # Check that stored velocities vn, vs, vL are consistent.
    @testset "Check stored velocities" begin
        local (; vs, vn, vL, fs,) = iter
        local (; α, α′,) = forcing
        for i ∈ eachindex(fs)
            vL_test = similar(vs[i])
            for j ∈ eachindex(vL_test)
                v⃗ₛ = vs[i][j]
                v⃗ₙ = vn[i][j]
                v⃗ₙₛ = v⃗ₙ - v⃗ₛ
                s⃗′ = fs[i][j, UnitTangent()]
                vL_test[j] = v⃗ₛ + α * s⃗′ × v⃗ₙₛ - α′ * s⃗′ × (s⃗′ × v⃗ₙₛ)
            end
            @test vL_test ≈ vL[i]
        end
    end

    @test impulse[begin][3] ≈ norm(impulse[begin])  # initially, the ring oriented in the Z direction
    @test impulse[end][3] ≈ norm(impulse[end])      # check that the orientation is the same at the end

    # In fact, vortex_impulse computes a normalised impulse, which for the case of a vortex
    # ring it's just its area A = π R².
    Rs = sqrt.(getindex.(impulse, 3) ./ T(π))
    @test T(2π) .* Rs ≈ length  # this verifies that the ring preserves its circular shape

    # We lost line length and energy
    @test length[end] < length[begin]
    @test energy[end] < energy[begin]

    plt = lineplot(
        time, length ./ length[begin]; title = "Ring with friction",
        name = "Length", xlabel = "Time", ylabel = "Relative change",
    )
    lineplot!(plt, time, energy ./ energy[begin]; name = "Energy")
    display(plt)

    iter
end

function test_ring_friction(::Type{T}) where {T}
    params = generate_biot_savart_parameters(T)
    S = define_curve(Ring(); translate = π, scale = π / 4)
    N = 48
    f = Filaments.init(S, ClosedFilament{T}, N, QuinticSplineMethod())

    # Create normal fluid forcing
    vn = @inferred UniformVectorField(zero(Vec3{T}))  # quiescent normal fluid
    α::T = 0.2    # Magnus force
    α′::T = 0.02  # drag force
    forcing = @inferred NormalFluidForcing(vn; α, α′)

    @testset "Static" test_ring_friction_static(f, params, forcing)
    @testset "Dynamic" test_ring_friction_dynamic(f, params, forcing)

    nothing
end

@testset "Ring friction ($T)" for T ∈ (Float32, Float64)
    test_ring_friction(T)
end
