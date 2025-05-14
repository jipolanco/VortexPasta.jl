using VortexPasta
using VortexPasta.Filaments
using VortexPasta.BiotSavart
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

function generate_biot_savart_parameters(::Type{T}) where {T}
    Γ = 1.0
    a = 1e-6
    Δ = 1/4
    L = 2π
    aspect = (1, 1, 2)
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
    )
end

function filaments_from_functions(funcs::NTuple{N, F}, args...) where {N, F <: Function}
    isempty(funcs) && return ()
    S = first(funcs)
    f = Filaments.init(S, args...)
    f_next = filaments_from_functions(Base.tail(funcs), args...)
    (f, f_next...)
end

function generate_filaments(N; Ls)
    p = PeriodicLine()  # unperturbed straight line
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

function simulate(prob::VortexFilamentProblem, forcing; dt_factor = 1.0, affect_t! = Returns(nothing))
    δ = minimum(node_distance, prob.fs)
    dt = BiotSavart.kelvin_wave_period(prob.p, δ) * dt_factor
    reconnect = ReconnectBasedOnDistance(δ / 2)
    refinement = RefineBasedOnSegmentLength(0.8 * δ)
    adaptivity = AdaptBasedOnSegmentLength(dt_factor) | AdaptBasedOnVelocity(δ / 2)
    iter = @inferred init(
        prob, RK4();
        forcing, affect_t!,
        dt, reconnect, refinement, adaptivity,
    )
    if affect_t! !== Returns(nothing)
        @test match(r"affect_t\!: Function", repr(iter)) !== nothing  # check that affect_t! appears in println(iter)
    end
    if forcing isa FourierBandForcing
        @test iter.vs !== iter.vL
    end
    E_init = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
    E_final = E_init
    ks_spec, Ek_init = Diagnostics.energy_spectrum(iter)
    while iter.t < prob.tspan[2]
        E = Diagnostics.kinetic_energy(iter; quad = GaussLegendre(3))
        E_final = E
        Nf = length(iter.fs)
        Np = sum(length, iter.fs; init = 0)
        # @show iter.t, E, Nf, Np
        step!(iter)
    end
    _, Ek_final = Diagnostics.energy_spectrum(iter)
    E_ratio = E_final / E_init
    spectra = (; ks = ks_spec, Ek_init, Ek_final,)
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
    @test data_bs.state.smoothed == true  # this means that the field needs to be unsmoothed (by reverting Gaussian filter)
    vs_hat = map(Array, data_bs.field)  # GPU -> CPU copy in case `field` is a tuple of GPU arrays
    ks_grid = data_bs.wavenumbers
    SyntheticFields.from_fourier_grid!(vs_band, vs_hat, ks_grid)  # copy coefficients of smoothed velocity field

    # Revert Gaussian filter
    α_ewald = cache_bs.params.α
    let (; cs, qs, Δks) = vs_band
        for i in eachindex(cs, qs)
            k⃗ = Δks .* qs[i]
            k² = sum(abs2, k⃗)
            cs[i] = cs[i] * exp(k² / (4 * α_ewald^2))  # deconvolve
        end
    end

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
    plt = @views lineplot(ks[2:end], Ek_init[2:end]; title, xscale = :log10, yscale = :log10, ylim = (1e-3, 1e-1), name = "Initial")
    @views lineplot!(plt, ks[2:end], Ek_final[2:end]; name = "Final")
    display(plt)
end

@testset "Forcing" begin
    # Compute self-induced velocity of vortex filaments
    N = 32
    p = generate_biot_savart_parameters(Float64)
    fs = generate_filaments(N; Ls = p.Ls)
    cache_bs = BiotSavart.init_cache(p, fs)
    vs = map(similar ∘ nodes, fs)
    velocity_on_nodes!(vs, cache_bs, fs)
    vs_self = map(copy, vs)
    α_ewald = p.α

    # Obtain coarse-grained velocity field in Fourier space
    data = BiotSavart.get_longrange_field_fourier(cache_bs)
    vs_hat = data.field
    state = data.state
    @test state.quantity == :velocity
    @test state.smoothed == true

    # Initialise random normal fluid velocity in Fourier space
    rng = StableRNG(42)
    Ls = p.Ls

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
        save_files = false
        plots = false
        Random.seed!(rng, 42)
        SyntheticFields.init_coefficients!(rng, vn, vn_rms)  # randomly set non-zero Fourier coefficients of the velocity field
        tmax = 1.0
        prob = VortexFilamentProblem(fs, (zero(tmax), tmax), p)
        @testset "NormalFluidForcing" begin
            forcing = @inferred NormalFluidForcing(vn; α, α′)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 1.2136651822544542
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
            # @show E_ratio  # = 1.080...
            @test 1.04 < E_ratio < 1.12  # energy increases more slowly than with steady forcing, which makes sense!
        end
        @testset "FourierBandForcing" begin
            forcing = @inferred FourierBandForcing(vn; α, α′)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 3.617499405142961
            @test 3.5 < E_ratio < 3.7
            save_files && filaments_to_vtkhdf("forcing_band.vtkhdf", iter)
            plots && plot_spectra(spectra; title = "FourierBandForcing")
        end
        @testset "FourierBandForcing (filtered vorticity)" begin
            forcing = @inferred FourierBandForcing(vn; α, α′, filtered_vorticity = true)
            (; iter, E_ratio, spectra) = simulate(prob, forcing)
            # @show E_ratio  # = 2.5623375372116857
            @test 2.5 < E_ratio < 2.6
            save_files && filaments_to_vtkhdf("forcing_band_filtered_vorticity.vtkhdf", iter)
            plots && plot_spectra(spectra; title = "FourierBandForcing (filtered vorticity)")
        end
    end

end
