using VortexPasta
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.BiotSavart: PseudoGPU  # for testing only
using VortexPasta.PredefinedCurves
using VortexPasta.Forcing
using VortexPasta.FilamentIO
using VortexPasta.SyntheticFields
using VortexPasta.Timestepping
using VortexPasta.Diagnostics
using LinearAlgebra: ×
using UnicodePlots
using Random
using StableRNGs
using Test

function generate_biot_savart_parameters(::Type{T}; aspect, kws...) where {T}
    Γ = 1.0
    a = 1e-6
    Δ = 1/4
    L = 2π
    Ls = aspect .* L
    β = 3.0
    rcut = L / 2
    α = β / rcut
    kmax = 2α * β
    M = ceil(Int, kmax * L / π)
    Ns = aspect .* M
    ParamsBiotSavart(
        T;
        Γ, α, a, Δ, rcut, Ls, Ns,
        backend_short = NaiveShortRangeBackend(),
        backend_long = NonuniformFFTsBackend(σ = T(1.5), m = HalfSupport(4)),
        quadrature = GaussLegendre(3),
        kws...,
    )
end

function filaments_from_functions(funcs::NTuple{N, F}, args...) where {N, F <: Function}
    isempty(funcs) && return ()
    S = first(funcs)
    f = Filaments.init(S, args...)
    f_next = filaments_from_functions(Base.tail(funcs), args...)
    (f, f_next...)
end

function generate_filaments(N; Ls, small_scales = false)
    r(t) = if small_scales
        0.01 * cispi(2t) + 0.01 * cispi(N/4 * t)
    else
        0.01 * cispi(2t)
    end
    p = PeriodicLine(; r)  # perturbed straight line
    T = Float64
    method = CubicSplineMethod()
    funcs = (
        define_curve(p; scale = Ls, translate = (1 / 4, 1 / 4, 1 / 2) .* Ls, orientation = +1),
        define_curve(p; scale = Ls, translate = (3 / 4, 1 / 4, 1 / 2) .* Ls, orientation = -1),
        define_curve(p; scale = Ls, translate = (1 / 4, 3 / 4, 1 / 2) .* Ls, orientation = -1),
        define_curve(p; scale = Ls, translate = (3 / 4, 3 / 4, 1 / 2) .* Ls, orientation = +1),
    )
    collect(filaments_from_functions(funcs, ClosedFilament{T}, N, method))
end

function simulate(
        prob::VortexFilamentProblem, forcing, dissipation = NoDissipation();
        dt_factor = 1.0, callback = Returns(nothing), affect_t! = Returns(nothing),
    )
    δ = minimum(node_distance, prob.fs)
    dt = BiotSavart.kelvin_wave_period(prob.p, δ) * dt_factor
    reconnect = ReconnectBasedOnDistance(δ / 2)
    refinement = RefineBasedOnSegmentLength(0.8 * δ)
    adaptivity = AdaptBasedOnSegmentLength(dt_factor) | AdaptBasedOnVelocity(δ / 2)
    iter = @inferred init(
        prob, RK4();
        forcing, dissipation, affect_t!, callback,
        dt, reconnect, refinement, adaptivity,
    )
    if affect_t! !== Returns(nothing)
        @test match(r"affect_t\!: Function", repr(iter)) !== nothing  # check that affect_t! appears in println(iter)
    end
    if dissipation isa SmallScaleDissipationBS
        @test match(r"dissipation: SmallScaleDissipationBS", repr(iter)) !== nothing  # check that dissipation parameters are printed
    end
    if forcing isa FourierBandForcing
        @test iter.vs !== iter.vL
    elseif forcing isa FourierBandForcingBS && dissipation isa NoDissipation
        @test iter.vs !== iter.vL
        # Check that vf is separately stored and that it contains the right thing.
        @test iter.vL ≈ iter.vs + iter.vf rtol=1e-16
        # @test iter.vL == iter.vs + iter.vf
    end
    E_init = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
    E_final = E_init
    ks_spec, Ek_init = Diagnostics.energy_spectrum(iter)
    _, Hk_init = Diagnostics.helicity_spectrum(iter)
    while iter.t < prob.tspan[2]
        E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
        E_final = E
        Nf = length(iter.fs)
        Np = sum(length, iter.fs; init = 0)
        # @show iter.t, E, Nf, Np
        step!(iter)
    end
    _, Ek_final = Diagnostics.energy_spectrum(iter)
    _, Hk_final = Diagnostics.helicity_spectrum(iter)
    E_ratio = E_final / E_init
    spectra = (; ks = ks_spec, Ek_init, Ek_final, Hk_init, Hk_final,)
    (; iter, E_ratio, spectra,)
end

# Check that the Fourier band forcing computes the right thing.
function check_fourier_band_forcing(forcing::FourierBandForcing, f::AbstractFilament, vL_computed, vs_self, cache_bs::BiotSavart.BiotSavartCache)
    vn_band = (forcing.vn)::FourierBandVectorField
    vs_band = (@inferred similar(vn_band))::FourierBandVectorField
    @test length(vs_band) == length(vn_band)  # same number of Fourier coefficients
    @test vs_band.qs == vn_band.qs    # same wavevectors
    @test vs_band.Δks == vn_band.Δks  # same wavevectors

    # Copy coefficients from vs_hat
    data_bs = @inferred BiotSavart.get_longrange_field_fourier(cache_bs)
    @test data_bs.state.quantity == :velocity
    @test data_bs.state.smoothing_scale == 0  # the field has not been smoothed
    vs_hat = map(Array, data_bs.field)  # GPU -> CPU copy in case `field` is a tuple of GPU arrays
    ks_grid = data_bs.wavenumbers
    SyntheticFields.from_fourier_grid!(vs_band, vs_hat, ks_grid)  # copy coefficients of smoothed velocity field

    # Compute v_{ns} in Fourier band
    vns_band = similar(vn_band)
    @. vns_band.cs = vn_band.cs - vs_band.cs

    # Estimate vL
    vL = similar(vL_computed)
    for i in eachindex(vL, f)
        s⃗ = f[i]
        s⃗′ = f[i, UnitTangent()]
        v⃗ₙₛ = vns_band(s⃗)
        vL[i] = vs_self[i] + forcing.α * s⃗′ × v⃗ₙₛ - forcing.α′ * s⃗′ × (s⃗′ × v⃗ₙₛ)
    end

    let a = collect(vL), b = collect(vL_computed)   # the `collect` is just in case these are PaddedVectors
        @test a ≈ b
    end

    nothing
end

function filaments_to_vtkhdf(filename, iter)
    Ls = iter.prob.p.Ls
    write_vtkhdf(filename, iter.fs; refinement = 4, periods = Ls) do io
        io["velocity_self"] = iter.vs
        io["velocity_total"] = iter.vL
    end
end

function plot_spectra(spectra; title)
    (; ks, Ek_init, Ek_final,) = spectra
    plt = @views lineplot(ks[2:end], Ek_init[2:end]; title, xscale = :log10, yscale = :log10, ylim = (1e-4, 1e-1), name = "Initial")
    @views lineplot!(plt, ks[2:end], Ek_final[2:end]; name = "Final")
    display(plt)
end

function get_energy_at_normalised_wavevector(cache::BiotSavart.BiotSavartCache, q⃗::NTuple{N}) where {N}
    Ls = cache.params.Ls
    Δks = @. 2 * (π / Ls)
    k⃗ = q⃗ .* Δks  # unnormalised wavevector
    k² = sum(abs2, k⃗)

    (; field, wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache)

    # Find wavevector
    wavenumbers::NTuple{N, AbstractVector}
    inds_pos = map(wavenumbers, k⃗) do ks, k
        findfirst(==(k), ks)
    end
    @assert all(!isnothing, inds_pos)

    # Unsmoothed velocity at wavevector k⃗
    @assert state.quantity == :velocity
    @assert state.smoothing_scale == 0
    v̂ = Vec3(map(u -> u[inds_pos...], field))

    # Energy at wavevectors k⃗ and -k⃗ (assuming Hermitian symmetry) -> hence no division by 2!
    Ek = sum(abs2, v̂)

    Ek
end

@testset "Forcing" begin
    save_files = false
    plots = false

    # Compute self-induced velocity of vortex filaments
    N = 32
    aspect = (1, 1, 2)
    params = generate_biot_savart_parameters(Float64; aspect)
    fs = generate_filaments(N; Ls = params.Ls)
    cache_bs = BiotSavart.init_cache(params, fs)
    vs = map(similar ∘ nodes, fs)
    velocity_on_nodes!(vs, cache_bs, fs)
    vs_self = map(copy, vs)
    α_ewald = params.α

    # Obtain velocity field in Fourier space
    data = BiotSavart.get_longrange_field_fourier(cache_bs)
    vs_hat = data.field
    state = data.state
    @test state.quantity == :velocity
    @test state.smoothing_scale == 0

    # Initialise random normal fluid velocity in Fourier space
    rng = StableRNG(42)
    Ls = params.Ls

    @testset "Static" begin
        # For testing purposes, we set a kmax which is larger than the kmax of the long-range grid.
        vn = @inferred FourierBandVectorField(undef, Ls; kmin = 0.1, kmax = 10.5)
        vn_rms = 1.0
        SyntheticFields.init_coefficients!(rng, vn, vn_rms)  # randomly set non-zero Fourier coefficients of the velocity field
        @testset "NormalFluidForcing" begin
            forcing = @inferred NormalFluidForcing(vn; α = 0.8, α′ = 0.2)
            cache = @inferred Forcing.init_cache(forcing, cache_bs)
            Forcing.update_cache!(cache, forcing, cache_bs)  # doesn't do anything for NormalFluidForcing
            @test startswith("NormalFluidForcing")(repr(forcing))
            for i in eachindex(fs, vs)
                copyto!(vs[i], vs_self[i])  # self-induced velocity (before forcing)
                Forcing.apply!(forcing, cache, vs[i], fs[i])
            end
        end
        @testset "FourierBandForcing" begin
            forcing = @inferred FourierBandForcing(vn; α = 0.8, α′ = 0.2)
            cache = @inferred Forcing.init_cache(forcing, cache_bs)
            Forcing.update_cache!(cache, forcing, cache_bs)
            @test startswith("FourierBandForcing")(repr(forcing))
            for i in eachindex(fs, vs)
                copyto!(vs[i], vs_self[i])  # self-induced velocity (before forcing)
                Forcing.apply!(forcing, cache, vs[i], fs[i])
                check_fourier_band_forcing(forcing, fs[i], vs[i], vs_self[i], cache_bs)
            end
        end
    end

    # We perform a short simulation and simply check that energy augmented.
    @testset "Simulation" begin
        vn_rms = 0.5
        α = 5.0
        α′ = 0.0
        vn = @inferred FourierBandVectorField(undef, Ls; kmin = 0.1, kmax = 1.5)
        Random.seed!(rng, 42)
        SyntheticFields.init_coefficients!(rng, vn, vn_rms)  # randomly set non-zero Fourier coefficients of the velocity field
        tmax = 1.0
        prob = VortexFilamentProblem(fs, (zero(tmax), tmax), params)
        @testset "NormalFluidForcing" begin
            forcing = @inferred NormalFluidForcing(vn; α, α′)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 1.2096699608685284
            @test 1.15 < E_ratio < 1.25
            save_files && filaments_to_vtkhdf("forcing_normal.vtkhdf", iter)
            plots && plot_spectra(spectra; title = "NormalFluidForcing")
        end
        @testset "NormalFluidForcing (time varying)" begin
            forcing = @inferred NormalFluidForcing(vn; α, α′)
            cs₀ = copy(vn.cs)  # initial Fourier coefficients of the forcing
            function affect_t!(iter, t)
                # Modify forcing with a chosen frequency
                ω = 1 * 2π / tmax
                @. iter.forcing.vn.cs = cs₀ * cis(ω * t)  # multiply initial coefficients by exp(im * ω * t)
            end
            (; iter, E_ratio, spectra) = simulate(prob, forcing; affect_t!)
            # @show E_ratio  # = 1.0795...
            @test 1.04 < E_ratio < 1.12  # energy increases more slowly than with steady forcing, which makes sense!
        end
        @testset "FourierBandForcing" begin
            forcing = @inferred FourierBandForcing(vn; α, α′)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 3.4612641897024696
            @test 3.4 < E_ratio < 3.8
            save_files && filaments_to_vtkhdf("forcing_band.vtkhdf", iter)
            plots && plot_spectra(spectra; title = "FourierBandForcing")
        end
        @testset "FourierBandForcing (filtered vorticity)" begin
            forcing = @inferred FourierBandForcing(vn; α, α′, filtered_vorticity = true)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 2.5823657991737927
            @test 2.45 < E_ratio < 2.68
            save_files && filaments_to_vtkhdf("forcing_band_filtered_vorticity.vtkhdf", iter)
            plots && plot_spectra(spectra; title = "FourierBandForcing (filtered vorticity)")
        end
        @testset "FourierBandForcingBS (constant α)" begin
            forcing = @inferred FourierBandForcingBS(; α = 100.0, kmin = 0.5, kmax = 2.5)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 1.1410585242365534
            @test 1.10 < E_ratio < 1.20
        end
        @testset "FourierBandForcingBS (constant α and α′)" begin
            forcing = @inferred FourierBandForcingBS(; α = 100.0, α′ = 1.0, kmin = 0.5, kmax = 2.5)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 1.1517750363770172
            @test 1.10 < E_ratio < 1.20
        end
        @testset "FourierBandForcingBS (constant ε_target, α′ = $α′)" for α′ in (0.0, 1.0)
            forcing = @inferred FourierBandForcingBS(; ε_target = 2e-2, α′, kmin = 1.5, kmax = 2.5)
            times = Float64[]
            energy_k = Float64[]  # energy at forced wavevectors
            energy = Float64[]
            ε_inj_time = Float64[]  # we currently don't use this
            function callback(iter)
                iter.nstep == 0 && empty!.((times, energy, energy_k))
                Ek = sum(iter.forcing_cache.qs) do q⃗  # sum energies of all forced wavevectors
                    get_energy_at_normalised_wavevector(iter.cache_bs, q⃗)
                end
                E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
                ε = Diagnostics.energy_injection_rate(iter; quad = GaussLegendre(3))
                push!(times, iter.t)
                push!(energy_k, Ek)
                push!(energy, E)
                push!(ε_inj_time, ε)
                nothing
            end
            (; iter, E_ratio, spectra) = simulate(prob, forcing; callback, dt_factor = 0.5)
            # @show α′, E_ratio  # = (0.0, 2.106383252320079) / (1.0, 2.100681284433628)
            @test 2.0 < E_ratio < 2.2
            let
                ks, fluxes = @inferred Diagnostics.energy_flux(iter, 8; quad = GaussLegendre(3))
                @test length(ks) ≤ 8
                @test length(fluxes) == 2
                @test hasproperty(fluxes, :vs)
                @test hasproperty(fluxes, :vf)
                @test fluxes.vf[end] > 0  # positive energy injection
            end
            # In the plot one should see that energy increases nearly linearly with the
            # imposed ε, until it starts to saturate near the end.
            a, b = 0.1, 0.4
            plots && let plt = lineplot(times, energy_k; xlim = (0, 1), ylim = (0, 0.02))
                lineplot!(plt, times, times * forcing.ε_target)
                vline!(plt, a)
                vline!(plt, b)
                display(plt)
            end
            let a = searchsortedlast(times, a), b = searchsortedlast(times, b)
                local ε_inj = (energy_k[b] - energy_k[a]) / (times[b] - times[a])  # mean energy injection rate at forced wavevectors
                # @show ε_inj
                @test ε_inj ≈ forcing.ε_target rtol=1e-2
            end
        end
        @testset "FourierBandForcingBS (single k⃗)" begin
            qs = [(1, 2, 3)]  # force a single wavevector
            forcing = @inferred FourierBandForcingBS(; ε_target = 1e-2, qs)
            # @test repr(forcing)
            times = Float64[]
            energy_k = Float64[]
            energy = Float64[]
            ε_inj_time = Float64[]  # we currently don't use this
            function callback(iter)
                iter.nstep == 0 && empty!.((times, energy, energy_k, ε_inj_time))
                Ek = get_energy_at_normalised_wavevector(iter.cache_bs, qs[1])
                E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
                ε = Diagnostics.energy_injection_rate(iter; quad = GaussLegendre(3))
                push!(times, iter.t)
                push!(energy_k, Ek)
                push!(energy, E)
                push!(ε_inj_time, ε)
                nothing
            end
            (; iter, E_ratio, spectra) = simulate(prob, forcing; callback)
            # In the following plot, one sees that energy at the active wavevector
            # linearly increases with the wanted ε_target at the beginning. Later it
            # saturates and tends to decrease, perhaps due to energy transfers to other
            # wavevectors.
            plots && let plt = lineplot(times, energy_k)
                lineplot!(plt, times, times * forcing.ε_target)
                display(plt)
            end
            save_files && open("energy_k.dat", "w") do io
                for i in eachindex(times, energy_k)
                    println(io, times[i], '\t', energy_k[i])
                end
            end
            # Check linear evolution at beginning of simulation, with the expected ε.
            let a = eachindex(times)[begin + 1], b = eachindex(times)[end ÷ 3]
                local ε_inj = (energy_k[b] - energy_k[a]) / (times[b] - times[a])  # energy injection rate at wavenumber k⃗
                # @show ε_inj
                @test ε_inj ≈ forcing.ε_target rtol=1e-3  # the ε_target really represents the energy injection rate at this wavevector!
            end
        end
    end

    @testset "Dissipation" begin
        tmax = 1.0
        fs_diss = generate_filaments(N; Ls, small_scales = true)
        prob_diss = VortexFilamentProblem(fs_diss, (zero(tmax), tmax), params)
        kdiss = 3  # dissipate at relatively small k for testing purposes

        times = Float64[]
        energy = Float64[]
        ε_inj_time = Float64[]
        ε_diss_time = Float64[]

        function callback(iter)
            iter.nstep == 0 && empty!.((times, energy, ε_inj_time, ε_diss_time))
            E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
            if iter.forcing isa FourierBandForcingBS
                ε_inj = Diagnostics.energy_injection_rate(iter, iter.vf; quad = GaussLegendre(3))
            else
                ε_inj = 0.0  # iter.vf is not defined
            end
            ε_diss = -Diagnostics.energy_injection_rate(iter, iter.vdiss; quad = GaussLegendre(3))
            push!(times, iter.t)
            push!(energy, E)
            push!(ε_inj_time, ε_inj)
            push!(ε_diss_time, ε_diss)
            nothing
        end

        @testset "DissipationBS: Constant α" begin
            dissipation = @inferred DissipationBS(; α = 0.1)
            (; iter, E_ratio, spectra) = simulate(prob_diss, NoForcing(), dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss
            @test all(>(0), ε_diss_time)  # energy dissipation is positive
            @test all(<(0), diff(energy) ./ diff(times))  # energy is actually lost
            plots && plot_spectra(spectra; title = "DissipationBS (constant α)")
        end

        @testset "DissipationBS: Constant ε_target" begin
            dissipation = @inferred DissipationBS(; ε_target = 1e-4)
            (; iter, E_ratio, spectra) = simulate(prob_diss, NoForcing(), dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss
            @test all(>(0), ε_diss_time)  # energy dissipation is positive
            @test all(<(0), diff(energy) ./ diff(times))  # energy is actually lost
            # Chosen dissipation rate corresponds almost exactly to actual dissipation rate.
            @test all(ε_d -> isapprox(ε_d, dissipation.ε_target; rtol = 2e-3), ε_diss_time)
            plots && plot_spectra(spectra; title = "DissipationBS (constant ε_target)")
            let
                ks, fluxes = @inferred Diagnostics.energy_flux(iter, 8; quad = GaussLegendre(3))
                @test length(ks) ≤ 8
                @test length(fluxes) == 2
                @test hasproperty(fluxes, :vs)
                @test hasproperty(fluxes, :vdiss)
            end
        end

        @testset "SmallScaleDissipationBS: Constant α" begin
            dissipation = @inferred SmallScaleDissipationBS(; α = 0.1, kdiss)
            (; iter, E_ratio, spectra) = simulate(prob_diss, NoForcing(), dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss
            @test all(>(0), ε_diss_time)  # energy dissipation is positive
            @test all(<(0), diff(energy) ./ diff(times))  # energy is actually lost
            plots && plot_spectra(spectra; title = "SmallScaleDissipationBS (constant α)")
        end

        @testset "SmallScaleDissipationBS: Constant ε_target" begin
            dissipation = @inferred SmallScaleDissipationBS(; ε_target = 1e-4, kdiss)
            (; iter, E_ratio, spectra) = simulate(prob_diss, NoForcing(), dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss
            @test all(>(0), ε_diss_time)  # energy dissipation is positive
            @test all(<(0), diff(energy) ./ diff(times))  # energy is actually lost
            # Chosen dissipation rate corresponds quite well to actual dissipation rate (differences are up to 5%).
            @test all(ε_d -> isapprox(ε_d, dissipation.ε_target; rtol = 0.05), ε_diss_time)
            plots && plot_spectra(spectra; title = "SmallScaleDissipationBS (constant ε_target)")
        end

        @testset "Forcing + dissipation (DissipationBS)" begin
            dissipation = @inferred DissipationBS(; ε_target = 1e-4)
            forcing = @inferred FourierBandForcingBS(; ε_target = 2e-4, kmin = 0.1, kmax = 2.5)
            (; iter, E_ratio, spectra) = simulate(prob_diss, forcing, dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss + iter.vf
            @test all(ε_d -> isapprox(ε_d, dissipation.ε_target; rtol = 0.05), ε_diss_time)
            let
                ks, fluxes = @inferred Diagnostics.energy_flux(iter, 8; quad = GaussLegendre(3))
                @test length(ks) ≤ 8
                @test length(fluxes) == 3
                @test hasproperty(fluxes, :vs)
                @test hasproperty(fluxes, :vf)
                @test hasproperty(fluxes, :vdiss)
            end
            plots && plot_spectra(spectra; title = "Forcing + dissipation")
        end

        @testset "Forcing + dissipation (SmallScaleDissipationBS)" begin
            dissipation = @inferred SmallScaleDissipationBS(; ε_target = 1e-4, kdiss)
            forcing = @inferred FourierBandForcingBS(; ε_target = 2e-4, kmin = 0.1, kmax = 2.5)
            (; iter, E_ratio, spectra) = simulate(prob_diss, forcing, dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss + iter.vf
            @test all(ε_d -> isapprox(ε_d, dissipation.ε_target; rtol = 0.05), ε_diss_time)
            plots && plot_spectra(spectra; title = "Forcing + dissipation")
        end

        @testset "Forcing + dissipation (PseudoGPU, SmallScaleDissipationBS)" begin
            backend_long = NonuniformFFTsBackend(PseudoGPU(); σ = 1.5, m = HalfSupport(4))
            params_gpu = generate_biot_savart_parameters(Float64; aspect, backend_long)
            prob_gpu = VortexFilamentProblem(fs_diss, (zero(tmax), tmax), params_gpu)
            dissipation = @inferred SmallScaleDissipationBS(; ε_target = 1e-4, kdiss)
            forcing = @inferred FourierBandForcingBS(; ε_target = 2e-4, kmin = 0.1, kmax = 2.5)
            (; iter, E_ratio, spectra) = simulate(prob_gpu, forcing, dissipation; callback)
            @test iter.vL ≈ iter.vs + iter.vdiss + iter.vf
            @test all(ε_d -> isapprox(ε_d, dissipation.ε_target; rtol = 0.05), ε_diss_time)
            plots && plot_spectra(spectra; title = "Forcing + dissipation (PseudoGPU)")
        end

        @testset "SmallScaleDissipationBS: Error if kdiss is too large" begin
            dissipation = @inferred SmallScaleDissipationBS(; ε_target = 0.1, kdiss = 10000)
            @test_throws "dissipative wavenumber (kdiss = $(dissipation.kdiss)) is too large" simulate(prob_diss, NoForcing(), dissipation; callback)
        end
    end
end
