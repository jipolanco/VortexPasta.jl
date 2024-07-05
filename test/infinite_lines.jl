# This tests infinite (non-closed) filaments with different discretisation methods.

using Test
using StaticArrays
using Statistics: mean, std
using LinearAlgebra: norm, cholesky
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.FilamentIO
using VortexPasta.BiotSavart
using VortexPasta.BiotSavart: LongRangeCache
using VortexPasta.Diagnostics
using UnicodePlots: lineplot, lineplot!

# Initialise nearly straight vortex line with helicoidal perturbation.
function init_vortex_line(; x, y, Lz = 2π, sign, A = 0.0, k::Int = 1)
    tlims = (-0.5, 0.5)
    S(t) = SVector(
        x + A * cospi(2 * k * t),
        y + A * sinpi(2 * k * t),
        (0.5 + sign * t) * Lz,
    )
    offset = setindex(zero(S(0.0)), sign * Lz, 3)
    (; x, y, Lz, sign, tlims, S, offset,)
end

function estimate_power_law_exponent(xs_in, ys_in)
    xs = log.(xs_in)
    Y = log.(ys_in)
    Base.require_one_based_indexing(xs)
    # Fit line using ordinary least squares
    X = hcat(fill!(similar(xs), 1), xs)
    βs = cholesky(X'X) \ (X'Y)
    (exp(βs[1]), βs[2])  # (prefactor, exponent)
end

# Check that I/O preserves the end-to-end offset.
function test_infinite_line_io(fs, vs; periods)
    with_periods = any(!isnothing, periods)
    method = Filaments.discretisation_method(first(fs))
    T = eltype(eltype(eltype(fs)))
    @assert T <: AbstractFloat
    suffix = with_periods ? "_periods" : ""
    fname = "infinite$(suffix).vtkhdf"
    FilamentIO.write_vtkhdf(fname, fs; periods, refinement = 2) do io
        io["velocity"] = vs
    end
    function check_velocity(io)
        vs_read = @inferred read(io, "velocity", FilamentIO.PointData())
        for (u, v) ∈ zip(vs, vs_read)
            # Note that u == v can fail when these are PaddedArrays since their ghost cells
            # generally don't match (and they don't have to). So we ignore possible ghost
            # cells in the comparison.
            @test @views u[:] == v[:]
        end
    end
    fs_read = @inferred FilamentIO.read_vtkhdf(check_velocity, fname, T, method)
    @test fs_read == fs
    for (f, g) ∈ zip(fs, fs_read)
        @test end_to_end_offset(f) == end_to_end_offset(g)
    end
    nothing
end

# Compute an extended energy spectrum with a chosen number of Fourier modes.
function extended_energy_spectrum(cache_bs::BiotSavart.BiotSavartCache, fs, Ns::Dims{3})
    cache_in = cache_bs.longrange
    cache = similar(cache_in, Ns)
    BiotSavart.add_point_charges!(cache, fs)
    BiotSavart.compute_vorticity_fourier!(cache)
    Diagnostics.energy_spectrum(cache)  # computes energy spectrum from vorticity in Fourier space
end

# Compute kinetic energy included in truncated Fourier coefficients of the vorticity.
# This is the example in the docs of BiotSavart.compute_on_nodes!.
function truncated_kinetic_energy_from_vorticity(cache::LongRangeCache)
    (; wavenumbers, uhat, ewald_prefactor,) = cache.common
    with_hermitian_symmetry = wavenumbers[1][end] > 0  # this depends on the long-range backend
    γ² = ewald_prefactor^2  # = (Γ/V)^2 [prefactor not included in the vorticity]
    E = 0.0
    for I ∈ CartesianIndices(uhat)
        k⃗ = map(getindex, wavenumbers, Tuple(I))
        kx = k⃗[1]
        factor = (!with_hermitian_symmetry || kx == 0) ? 0.5 : 1.0
        k² = sum(abs2, k⃗)
        if !iszero(k²)
            ω⃗ = uhat[I]  # Fourier coefficient of the vorticity
            E += γ² * factor * sum(abs2, ω⃗) / k²
        end
    end
    E
end

function test_infinite_lines(method)
    L = 4π
    lines = let
        A = 0.08 * L / 2π
        k = 2
        h = L/4
        [
            init_vortex_line(; x = 1h, y = 1h, Lz = L, sign = +1, A, k,),
            init_vortex_line(; x = 1h, y = 3h, Lz = L, sign = -1, A, k,),
            init_vortex_line(; x = 3h, y = 1h, Lz = L, sign = -1, A, k,),
            init_vortex_line(; x = 3h, y = 3h, Lz = L, sign = +1, A, k,),
        ]
    end

    N = 32
    filaments = map(lines) do line
        (; tlims, S, offset,) = line
        ζs = range(tlims...; length = N + 1)[1:N]
        @inferred Filaments.init(ClosedFilament, S.(ζs), method; offset)
    end

    @testset "Filaments" begin
        f = first(filaments)
        update_coefficients!(f)
        Xoffset = end_to_end_offset(f)
        @test Xoffset[1] == Xoffset[2] == 0
        @test abs(Xoffset[3]) == L
        # Xoffset = Vec3(0, 0, 2π)
        ta, tb = knotlims(f)
        T = tb - ta
        t₀ = 0.1  # arbitrary location
        @test f(t₀) + Xoffset ≈ f(t₀ + T)
        @test f(t₀, Derivative(1)) ≈ f(t₀ + T, Derivative(1))
        @test f(t₀, Derivative(2)) ≈ f(t₀ + T, Derivative(2))
    end

    @testset "Evaluation on knots" begin
        f = first(filaments)
        update_coefficients!(f)
        disc = Filaments.discretisation_method(f)
        for i ∈ (2, 5, lastindex(f))
            @test f[i] ≈ f(i, 0.0)
            if Filaments.continuity(disc) ≥ 1
                @test f[i, Derivative(1)] ≈ f(i, 0.0, Derivative(1))
            end
            if Filaments.continuity(disc) ≥ 2
                @test f[i, Derivative(2)] ≈ f(i, 0.0, Derivative(2))
            end
        end
    end

    @testset "Change offset" begin
        f = first(filaments)
        fc = copy(f)
        off = Vec3(1, 1, 1)
        g = @inferred Filaments.change_offset(fc, off)
        for (u, v) ∈ zip(Filaments.allvectors.((fc, g))...)
            @test u === v  # we're reusing the same arrays (for nodes, coefficients, etc...)
        end
        update_coefficients!(g)
        @test end_to_end_offset(g) == off
        ta, tb = knotlims(g)
        T = tb - ta
        t₀ = 0.1  # arbitrary location
        @test g(t₀) + off ≈ g(t₀ + T)
    end

    Ls = (L, L, L)
    Ns = (1, 1, 1) .* 64
    kmax = (Ns[1] ÷ 2) * 2π / L
    α = kmax / 4
    rcut = 4 / α

    params = ParamsBiotSavart(;
        Γ = 2.0,
        a = 1e-6,
        Δ = 1/4,
        Ls, Ns, rcut, α,
        longrange_truncate_spherical = true,  # just for testing this option, shouldn't change much...
    )
    cache = @inferred BiotSavart.init_cache(params, filaments)
    vs = map(f -> similar(nodes(f)), filaments)
    ψs = map(f -> similar(nodes(f)), filaments)

    E_from_vorticity = Ref(0.0)  # large-scale kinetic energy from vorticity field

    function callback_vorticity(cache)
        E = truncated_kinetic_energy_from_vorticity(cache)
        E_from_vorticity[] = E
        nothing
    end

    compute_on_nodes!(
        (velocity = vs, streamfunction = ψs), cache, filaments;
        callback_vorticity,
    )

    @testset "Velocities" begin
        vnorms = norm.(first(vs))
        vmean = mean(vnorms)
        vstd = std(vnorms)
        # The velocity std should be approximately zero.
        # This condition fails for continuity == 0 (HermiteInterpolation(0), i.e. linear segments).
        continuity = Filaments.continuity(Filaments.interpolation_method(first(filaments)))
        if continuity ≥ 1
            @test vstd / abs(vmean) < 1e-3
        end
    end

    E = Diagnostics.kinetic_energy(filaments, ψs, params.Γ, params.Ls)

    @testset "Energy spectrum" begin
        ks, Ek = Diagnostics.energy_spectrum(cache)
        local (; Γ, Ls,) = params
        local Lvort = Diagnostics.filament_length(filaments; quad = GaussLegendre(4))
        local Cspec = Γ^2 * Lvort / (4π * prod(Ls))  # analytical prefactor of energy spectrum at large k
        Ns_ext = @. (Ns ÷ 2) * 3  # compute spectrum with 1.5× resolution
        ks_ext, Ek_ext = extended_energy_spectrum(cache, filaments, Ns_ext)
        @views plt = lineplot(
            ks[2:end], Ek[2:end];
            xscale = log10, yscale = log10,
            title = "Infinite lines", xlabel = "k", ylabel = "E(k)",
            name = "Spectrum", xlim = (0.8 * 2π/L, max(Ns_ext...) * 0.6),
            ylim = (1e-4, 4e-2),
        )
        @views lineplot!(plt, ks_ext[2:end], Ek_ext[2:end]; name = "Spectrum (extended)")
        @views let ks = ks[4:end]
            local ys = Cspec ./ ks
            lineplot!(plt, ks, ys; name = "Analytical (~k^{-1})")
        end
        @views let inds = (lastindex(ks) ÷ 4):(lastindex(ks) - 2)
            C, α = estimate_power_law_exponent(ks[inds], Ek[inds])
            @test isapprox(α, -1; rtol = 1e-2)     # approximately k^{-1}
            @test isapprox(C, Cspec; rtol = 0.04)  # the fitted prefactor is close to the analytical one
            lineplot!(plt, ks[inds], @.(0.5 * C * ks[inds]^α); name = "Fit (shifted)")
        end
        println(plt)

        # The extended spectrum should contain more energy than the original one, but less
        # energy than the total energy of the system (since in both cases we're discarding
        # energy at very small scales).
        Δk = ks[2] - ks[1]
        E_spec = sum(Ek) * Δk
        E_spec_ext = sum(Ek_ext) * Δk

        # @show E E_spec_ext E_spec E_from_vorticity[] E_spec_ext/E E_spec/E
        # @show (E_spec - E_from_vorticity[]) / E_spec
        @test isapprox(E_spec, E_from_vorticity[]; rtol = 1e-6)
        @test E_spec < E_spec_ext < E

        # Check that we're in the right order of magnitude (actual limits will depend on
        # parameters...).
        @test 0.2 < E_spec/E < 0.3
        @test 0.2 < E_spec_ext/E < 0.3
    end

    @testset "FilamentIO" begin
        test_infinite_line_io(filaments, vs; periods = (nothing, nothing, nothing))
        test_infinite_line_io(filaments, vs; periods = Ls)
    end

    nothing
end

@testset "Infinite lines" begin
    methods = (
        "FiniteDiff(2) / Hermite(2)" => FiniteDiffMethod(2, HermiteInterpolation(2)),
        "FiniteDiff(2) / Hermite(1)" => FiniteDiffMethod(2, HermiteInterpolation(1)),
        "FiniteDiff(2) / Hermite(0)" => FiniteDiffMethod(2, HermiteInterpolation(0)),
        "CubicSpline" => CubicSplineMethod(),
    )
    @testset "$name" for (name, method) ∈ methods
        test_infinite_lines(method)
    end
end
