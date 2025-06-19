using Test
using VortexPasta
using VortexPasta.Filaments
using VortexPasta.PredefinedCurves
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Diagnostics: Diagnostics
using UnicodePlots: UnicodePlots, lineplot, lineplot!
using Rotations: Rotations, QuatRotation
using Random
using StableRNGs

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

function generate_biot_savart_parameters(::Type{T}, L::T) where {T}
    Γ = one(T)
    a = T(1e-6)
    Δ = T(1/4)
    Ngrid = 16
    backend_short = CellListsBackend(2)
    Ls = (L, L, L)
    Ns = (Ngrid, Ngrid, Ngrid)
    kmax::T = (Ngrid ÷ 2) * 2π / L
    β = 3
    α = kmax / 2β
    rcut::T = β / α
    backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4))
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short,
        backend_long,
        quadrature = GaussLegendre(3),
    )
end

function generate_random_ring(rng::AbstractRNG, ::Type{T}, Np, method; Ls, R) where {T}
    p = Ring()
    scale = R
    rotate = rand(rng, QuatRotation{T})
    translate = rand(rng, Vec3{T}) .* Vec3(Ls)
    S = define_curve(p; scale, rotate, translate)
    Filaments.init(S, ClosedFilament{T}, Np, method)
end

function test_injection(::Type{T}, Np, method) where {T}
    L = T(2) * π
    Ls = (L, L, L)
    R = L/2

    # Start with 4 random rings
    rng = StableRNG(42)
    fs = [generate_random_ring(rng, T, Np, method; Ls, R) for _ ∈ 1:4]
    d_min = minimum_node_distance(fs) * T(0.8)

    params = generate_biot_savart_parameters(T, L)
    (; Γ,) = params
    τ = R^2 / Γ / 2
    tspan = (zero(τ), τ)

    prob = VortexFilamentProblem(fs, tspan, params)

    dt = BiotSavart.kelvin_wave_period(params, d_min)
    scheme = RK4()
    reconnect = ReconnectBasedOnDistance(d_min)
    refinement = RefineBasedOnSegmentLength(d_min, 2 * d_min)
    adaptivity = AdaptBasedOnSegmentLength(0.8) | AdaptBasedOnVelocity(1.0)

    time = T[]
    length = T[]
    energy = T[]

    function callback(iter)
        local (; nstep, t,) = iter
        local quad = GaussLegendre(3)
        local E = Diagnostics.kinetic_energy(iter; quad)
        local L = Diagnostics.filament_length(iter; quad)
        # @show nstep, t, E, L
        push!(time, t)
        push!(length, L)
        push!(energy, E)
        nothing
    end

    inject_dt = τ / 12  # injection period
    inject_t_next = Ref(tspan[1])  # next time at which to inject (initially at the initial time, tspan[1] = 0)

    function affect!(iter)
        local (; t,) = iter
        if t ≥ inject_t_next[]
            local f = generate_random_ring(rng, T, Np, method; Ls, R)
            inject_filament!(iter, f)
            inject_t_next[] += inject_dt  # next time at which to inject
        end
        nothing
    end

    iter = @inferred init(
        prob, scheme;
        dt,
        adaptivity,
        refinement,
        reconnect,
        callback, affect!,
    )

    solve!(iter)

    if VERBOSE
        plt = lineplot(
            time, length ./ length[begin]; title = "Filament injection",
            name = "Length", xlabel = "Time", ylabel = "Relative change",
        )
        lineplot!(plt, time, energy ./ energy[begin]; name = "Energy")
        display(plt)
    end

    L_relchange = length[end] / length[begin]
    VERBOSE && @show L_relchange

    # Note: results can vary from one run to another, especially when running with multiple
    # threads. This is probably influenced by parallelisation of reconnections?
    @test 1.8 < L_relchange < 2.3  # we're at roughly 2× the initial vortex length

    nothing
end

@testset "Filament injection" begin
    Np = 32
    test_injection(Float64, Np, CubicSplineMethod())
end
