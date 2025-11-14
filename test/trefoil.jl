using Test
using LinearAlgebra: norm, ⋅, ×
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot
using VortexPasta.Filaments
using VortexPasta.Filaments: GeometricQuantity, Vec3
using VortexPasta.BiotSavart
using VortexPasta.Diagnostics
using ForwardDiff: ForwardDiff
using FFTW: FFTW
using LaTeXStrings  # used for plots only (L"...")

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

function trefoil_function()
    R = π / 3
    define_curve(TrefoilKnot(); translate = R, scale = R)
end

function trefoil_derivative(::Derivative{1}, t)
    S = trefoil_function()
    ForwardDiff.derivative(S, t)
end

function trefoil_derivative(::Derivative{n}, t) where {n}
    ForwardDiff.derivative(u -> trefoil_derivative(Derivative(n - 1), u), t)
end

# Analytical trefoil torsion (here t ∈ [0, 1]).
function trefoil_quantity(t, ::TorsionScalar)
    s′ = trefoil_derivative(Derivative(1), t)
    s″ = trefoil_derivative(Derivative(2), t)
    s‴ = trefoil_derivative(Derivative(3), t)
    b⃗ = s′ × s″
    (b⃗ ⋅ s‴) / norm(b⃗)^2
end

function trefoil_quantity(t, ::CurvatureScalar)
    s′ = trefoil_derivative(Derivative(1), t)
    s″ = trefoil_derivative(Derivative(2), t)
    s′² = sum(abs2, s′)
    norm(s′² * s″ - (s″ ⋅ s′) * s′) / (s′² * s′²)
end

nderivs_required(::GeometricQuantity) = Val(2)
nderivs_required(::TorsionScalar) = Val(3)

function trefoil_quantity_error(N, method, quantity)
    f = init_trefoil_filament(
        N; method, nderivs = nderivs_required(quantity),
    )
    τ_max = 0.0
    err = 0.0
    err_max = 0.0
    for i ∈ eachindex(f)
        τ = f[i, quantity]
        u = (i - 1) / N  # original curve parametrisation
        τ_expected = trefoil_quantity(u, quantity)
        τ_max = max(τ_max, abs(τ_expected))
        err_i = (τ - τ_expected)^2
        err += err_i
        err_max = max(err_max, err_i)
    end
    err = sqrt(err / N) / τ_max  # normalise by the maximum torsion
    err_max = sqrt(err_max) / τ_max
    err, err_max
end

# This tests scalar quantities on discretisation points, comparing them with "analytical"
# expressions obtained via automatic differentiation (using ForwardDiff).
function test_trefoil_quantity(method, quantity)
    N = 96
    err, err_max = trefoil_quantity_error(N, method, quantity)
    # @show err, err_max
    if method isa QuinticSplineMethod
        # This method captures the torsion with a quite good precision.
        # Of course, the error decreases with the resolution `N` (with polynomial convergence).
        # In both cases, the error decreases as N⁻⁴.
        if quantity === TorsionScalar()
            @test err < 2e-4
            @test err_max < 3e-4
        elseif quantity === CurvatureScalar()
            @test err < 1e-4
            @test err_max < 2e-4
        end
    elseif method isa CubicSplineMethod
        # This error decreases as N⁻².
        @assert quantity === CurvatureScalar()
        @test err < 5e-3
        @test err_max < 1e-2
    elseif method isa FiniteDiffMethod
        # This error decreases as N⁻¹.
        @assert quantity === CurvatureScalar()
        @test err < 2e-2
        @test err_max < 3e-2
    elseif method isa FourierMethod
        # This is usually the case of FourierMethod, which captures the torsion with very
        # high precision! This is expected because the trefoil is defined by trigonometric
        # functions, so that very few Fourier modes are needed to perfectly describe the
        # trefoil. What is less expected is that this error tends to *increase* with `N`
        # (but always staying very small). This seems to be due to round-off error of the
        # FFTs. Indeed, Fourier coefficients that are supposed to be 0 are actually ~1e-16,
        # and are amplified when computing derivatives.
        if quantity === TorsionScalar()
            @test err < 1.1e-12
            @test err_max < 1e-11
        elseif quantity === CurvatureScalar()
            # Case of CurvatureScalar
            @test err < 1e-13
            @test err_max < 2e-13
        end
    end
    nothing
end

function test_convergence(method, quantity)
    Ns = 1 .<< (3:10)  # = 2^n
    errs_both = map(Ns) do N
        err, err_max = trefoil_quantity_error(N, method, quantity)
        err, err_max
    end
    errs = getindex.(errs_both, 1)
    errs_max = getindex.(errs_both, 2)
    Ns, errs, errs_max
end

function init_trefoil_filament(N::Int; method = CubicSplineMethod(), kws...)
    S = @inferred trefoil_function()
    ζs = range(0, 1; length = N + 1)[1:N]
    Xs = @inferred broadcast(S, ζs)  # same as S.(ζs)
    Filaments.init(ClosedFilament, Xs, method; kws...)
end

function compare_long_range(
        fs::AbstractVector{<:AbstractFilament}, backend;
        tol = 1e-8, params_kws...,
    )
    params_exact = @inferred ParamsBiotSavart(;
        params_kws...,
        backend_long = ExactSumBackend(),
    )
    params_default = @inferred ParamsBiotSavart(;
        params_kws...,
        backend_long = backend,
    )
    quad = params_exact.quad
    T = eltype(params_exact)

    # This is just to check that Base.show is implemented for ParamsBiotSavart.
    @test startswith(repr(params_exact), "ParamsBiotSavart{$T} with:\n")
    @test startswith(repr(params_default), "ParamsBiotSavart{$T} with:\n")

    cache_exact_base = @inferred(BiotSavart.init_cache(params_exact, fs))
    cache_default_base = @inferred(BiotSavart.init_cache(params_default, fs))

    cache_exact = cache_exact_base.longrange
    cache_default = cache_default_base.longrange

    @test BiotSavart.backend(cache_exact) isa ExactSumBackend
    @test BiotSavart.backend(cache_default) === backend

    # Compute induced velocity field in Fourier space
    foreach((cache_exact, cache_default)) do c
        BiotSavart.add_point_charges!(c.common.pointdata_d, fs, quad)
        BiotSavart.compute_vorticity_fourier!(c)
        BiotSavart.compute_velocity_fourier!(c)
    end

    # Compare velocities in Fourier space.
    ks_default = cache_default.common.wavenumbers_d
    ks_exact = cache_exact.common.wavenumbers_d
    inds_to_compare = ntuple(Val(3)) do i
        inds = eachindex(ks_exact[i])
        js = i == 1 ? inds[begin:end - 1] : inds
        @assert @views ks_default[i][js] == ks_exact[i][js]  # wavenumbers match in this index region
        js
    end |> CartesianIndices
    diffnorm_L2 = sum(inds_to_compare) do I
        uhat, vhat = cache_exact.common.uhat_d, cache_default.common.uhat_d
        du = uhat[I] - vhat[I]
        sum(abs2, du)  # = |u - v|^2
    end
    @test diffnorm_L2 < tol^2

    # Interpolate long-range velocity back to filament positions
    foreach((cache_exact, cache_default)) do c
        BiotSavart.set_interpolation_points!(c, fs)
        local callback = BiotSavart.get_ewald_interpolation_callback(c)  # this is to make sure we get the *long-range* (smoothed) velocity
        BiotSavart.interpolate_to_physical!(callback, c)
    end

    charges(cache) = cache.common.pointdata_d.charges  # get charges from cache
    max_rel_error_physical = maximum(zip(charges(cache_exact), charges(cache_default))) do (qexact, qdefault)
        norm(qexact - qdefault) / norm(qexact)
    end
    @test max_rel_error_physical < tol

    # Copy data to arrays.
    vs_exact = map(similar ∘ nodes, fs)
    vs_default = map(similar ∘ nodes, fs)
    BiotSavart.copy_long_range_output!(vs_exact, cache_exact)
    BiotSavart.copy_long_range_output!(vs_default, cache_default)

    # Compare velocities one filament at a time.
    function l2_difference(us, vs)
        T = eltype(eltype(eltype(us)))
        err = zero(T)
        n = 0
        for (u, v) ∈ zip(us, vs)
            for i ∈ eachindex(u, v)
                err += sum(abs2, u[i] - v[i]) / sum(abs2, v[i])
                n += 1
            end
        end
        sqrt(err / n)
    end

    err_velocity_long = l2_difference(vs_default, vs_exact)
    @test err_velocity_long < 2 * tol

    # Also compare output of full (short + long range) computation, and including
    # the streamfunction.
    ψs_exact = map(similar ∘ nodes, fs)
    ψs_default = map(similar ∘ nodes, fs)

    fields_exact = (;
        velocity = vs_exact,
        streamfunction = ψs_exact,
    )
    fields_default = (;
        velocity = vs_default,
        streamfunction = ψs_default,
    )

    BiotSavart.compute_on_nodes!(fields_default, cache_default_base, fs)
    BiotSavart.compute_on_nodes!(fields_exact, cache_exact_base, fs)

    err_velocity = l2_difference(vs_default, vs_exact)
    err_streamfunction = l2_difference(ψs_default, ψs_exact)
    @test err_velocity < tol
    @test err_streamfunction < tol

    nothing
end

function compare_short_range(fs::AbstractVector{<:AbstractFilament}; params_kws...)
    params_naive = @inferred ParamsBiotSavart(;
        params_kws...,
        backend_short = NaiveShortRangeBackend(),
    )
    params_cl = @inferred ParamsBiotSavart(;
        params_kws...,
        backend_short = CellListsBackend(2),
    )

    let params = params_naive
        rc_max = @inferred BiotSavart.max_cutoff_distance(params.shortrange.backend, params.Ls) 
        Lmin = min(params.Ls...)
        @test rc_max ≈ Lmin / 2
    end
    let params = params_cl
        rc_max = @inferred BiotSavart.max_cutoff_distance(params.shortrange.backend, params.Ls)
        Lmin = min(params.Ls...)
        @test Lmin / 3 < rc_max < Lmin / 2  # the actual value depends on the nsubdiv parameter (here nsubdiv = 2)
    end

    cache_naive = @inferred BiotSavart.init_cache(params_naive, fs)
    cache_cl = @inferred BiotSavart.init_cache(params_cl, fs)

    vs_naive = @inferred BiotSavart.velocity_on_nodes(cache_naive, fs; longrange = false)
    vs_cl = @inferred BiotSavart.velocity_on_nodes(cache_cl, fs; longrange = false)

    for (a, b) ∈ zip(vs_naive, vs_cl)
        # Note: since `a` and `b` are PaddedVectors, we need to exclude their ghost cells
        # (which can be different) from the comparison. For that reason we use a view onto
        # the "central" (non-ghost) data of both arrays.
        @test @views isapprox(a[:], b[:]; rtol = 1e-7)
    end

    nothing
end

function compute_filament_velocity_and_streamfunction(f::AbstractFilament; α, Ls, Ns, params_kws...)
    rcut = 4 * sqrt(2) / α
    @assert rcut < minimum(Ls) * 2 / 5  # cell lists requirement (with nsubdiv = 2)
    params = ParamsBiotSavart(; params_kws..., α, Ns, Ls, rcut)
    fs = [f]
    cache = init_cache(params, fs)
    vs = map(similar ∘ nodes, fs)
    fields = (
        velocity = vs,
        streamfunction = map(similar, vs),
    )
    compute_on_nodes!(fields, cache, fs)
    fields
end

# Check that the total induced velocity doesn't depend strongly on the Ewald parameter α.
# (In theory it shouldn't depend at all...)
function check_independence_on_ewald_parameter(f, αs; β = 3.5, quad = GaussLegendre(3), params_kws...)
    fields_all = map(αs) do α
        kmax = 2 * α * β
        Ns = ceil.(Int, (kmax / π) .* Ls) .+ 2
        compute_filament_velocity_and_streamfunction(
            f;
            α, Ns, rcut = β / α,
            backend_short = CellListsBackend(2),
            backend_long = NonuniformFFTsBackend(fftw_flags = FFTW.ESTIMATE),  # use FFTW.ESTIMATE to save some time
            quadrature = quad,
            params_kws...,
        )
    end
    fields_test = last(fields_all)
    maxdiffs_vel = map(fields_all) do fields
        maximum(zip(fields.velocity, fields_test.velocity)) do (a, b)
            norm(a - b) / norm(b)
        end
    end
    maxdiffs_stf = map(fields_all) do fields
        maximum(zip(fields.streamfunction, fields_test.streamfunction)) do (a, b)
            norm(a - b) / norm(b)
        end
    end
    # @show maxdiffs_vel maxdiffs_stf
    @test maximum(maxdiffs_vel) < 3e-4
    @test maximum(maxdiffs_stf) < 1e-5
    nothing
end

function test_helicity(f; quadrature = GaussLegendre(3), Ns, Ls, params_kws...)
    β = 3.5
    kmax = minimum(2π .* Ns ./ Ls)
    α = kmax / 2β
    rcut = β / α
    backend_short = CellListsBackend(2)
    params = ParamsBiotSavart(; α, Ns, Ls, rcut, quadrature, backend_short, params_kws...)
    fs = [f]
    cache = init_cache(params, fs)
    vs = map(similar ∘ nodes, fs)
    fields = (
        velocity = vs,
    )
    ks, Hk = Diagnostics.init_spectrum(cache)
    function callback_vorticity(cache::BiotSavart.LongRangeCache)
        Diagnostics.helicity_spectrum!(Hk, ks, cache)  # compute helicity spectrum directly from vorticity
        nothing
    end
    compute_on_nodes!(fields, cache, fs; callback_vorticity)
    Hk_from_vort = copy(Hk)
    fill!(Hk, 0)
    Diagnostics.helicity_spectrum!(Hk, ks, cache)  # now compute helicity spectrum from smoothed velocity
    @test Hk_from_vort ≈ Hk  # compare two different implementations
    H = Diagnostics.helicity(fs, vs, params; quad = quadrature)
    # open("Hk.dat", "w") do io
    #     for H in Hk
    #         println(io, H)
    #     end
    # end
    dk = ks[2] - ks[1]
    H_from_spectrum = sum(Hk) * dk
    # @show H H_from_spectrum
    @test 3H < H_from_spectrum < H < 0  # very rough approximation (spectrum decays very slowly; in fact I'm not sure it actually decays!)
    nothing
end

longrange_modified_bs_integrand_taylor(::Velocity, s⃗′, r⃗, α) =
    (s⃗′ × r⃗) * (4 * α^3 / (3 * sqrt(π)))

longrange_modified_bs_integrand_taylor(::Streamfunction, s⃗′, r⃗, α) =
    s⃗′ * (2 * α / sqrt(π))

function test_long_range_accuracy_near_zero(::Type{T}, quantity) where {T}
    component = BiotSavart.LongRange()
    s⃗′ = Vec3{T}(0.3, 0.4, -0.5)  # not important
    α = T(2.3)
    @testset "Near zero" begin
        r⃗ = Vec3{T}(0.0, 0.0, 0.0) .+ T(8) * sqrt(eps(1.0))
        r² = sum(abs2, r⃗)
        integrand = BiotSavart.modified_bs_integrand(quantity, component, s⃗′, r⃗, r², α)
        integrand_taylor = longrange_modified_bs_integrand_taylor(quantity, s⃗′, r⃗, α)
        @test isapprox(integrand, integrand_taylor; rtol = 1e-4)
    end
    @testset "At zero" begin
        r⃗ = Vec3{T}(0.0, 0.0, 0.0)
        r² = sum(abs2, r⃗)
        # This should call the internal BiotSavart.long_range_bs_integrand_at_zero function.
        integrand = BiotSavart.modified_bs_integrand(quantity, component, s⃗′, r⃗, r², α)
        integrand_taylor = longrange_modified_bs_integrand_taylor(quantity, s⃗′, r⃗, α)
        @test integrand == integrand_taylor
    end
    nothing
end

@testset "Trefoil" begin
    # Note: for testing purposes it's important to discretise the vortex with a *small*
    # number of discretisation points.
    f = @inferred init_trefoil_filament(30)
    Ls = (1.5π, 1.5π, 2π)  # Ly is small to test periodicity effects
    Ns = (3, 3, 4) .* 30
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    params_kws = (; Ls, Γ = 7.4, a = 1e-5,)
    β = 3.5  # accuracy coefficient
    α_default = 3.0
    kmax_default = 2 * α_default * β
    Ns_default = ceil.(Int, (kmax_default / π) .* Ls) .+ 2
    @testset "Long range" begin
        # Test NUFFT backends with default parameters
        tol = 1e-5  # allowed relative error
        @testset "NonuniformFFTsBackend" begin
            local β = 3.5
            local α = 3.0
            local kmax = 2 * α * β
            local Ns = ceil.(Int, (kmax / π) .* Ls) .+ 2
            compare_long_range([f], NonuniformFFTsBackend(fftw_flags = FFTW.ESTIMATE); tol, params_kws..., Ns, α, rcut = β / α)
        end
        @testset "$quantity near r = 0" for quantity ∈ (Velocity(), Streamfunction())
            test_long_range_accuracy_near_zero(Float64, quantity)
        end
    end
    @testset "Short range" begin
        compare_short_range([f]; params_kws..., Ns = Ns_default, α = α_default, rcut = β / α_default)
    end
    @testset "Dependence on α" begin
        αs = [kmax / 5, kmax / 6, kmax / 7, kmax / 8, kmax / 12, kmax / 16]
        quadratures = (GaussLegendre(3), NoQuadrature())
        @testset "$quad" for quad ∈ quadratures
            @testset "use_simd = $use_simd" for use_simd ∈ (true, false)
                @testset "avoid_explicit_erf = $avoid_explicit_erf" for avoid_explicit_erf ∈ (true, false)
                    check_independence_on_ewald_parameter(f, αs; quad, β, use_simd, avoid_explicit_erf, params_kws...)
                end
            end
        end
        @testset "FourierMethod()" begin
            f_fourier = @inferred init_trefoil_filament(32; method = FourierMethod())
            check_independence_on_ewald_parameter(f_fourier, αs; β, params_kws...)
        end
    end
    @testset "Helicity" begin
        test_helicity(f; Ns = Ns_default, α = α_default, rcut = β / α_default, params_kws...)
    end
    @testset "Curvature" begin
        methods = (FiniteDiffMethod(), CubicSplineMethod(), QuinticSplineMethod(), FourierMethod())
        @testset "$method" for method ∈ methods
            test_trefoil_quantity(method, CurvatureScalar())
        end
    end
    @testset "Torsion" begin
        methods = (QuinticSplineMethod(), FourierMethod())
        @testset "$method" for method ∈ methods
            test_trefoil_quantity(method, TorsionScalar())
        end
    end
end

##

if @isdefined(Makie)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    # wireframe!(ax, Rect(0, 0, 0, Ls...); color = :grey, linewidth = 0.5)
    plot!(ax, f; refinement = 8)
    fig
end

function plot_convergence()
    # Check convergence of discretisation
    fig = Figure(size = (1200, 600))
    ax_curv = Axis(fig[1, 1]; xscale = log10, yscale = log10, title = "Curvature", xlabel = L"N", ylabel = "Error")
    ax_tors = Axis(fig[1, 2];  xscale = log10, yscale = log10, title = "Torsion", xlabel = L"N")
    linkaxes!(ax_curv, ax_tors)
    methods_curv = (FiniteDiffMethod(), CubicSplineMethod(), QuinticSplineMethod(), FourierMethod())
    methods_tors = (QuinticSplineMethod(), FourierMethod())
    ls_curv = map(methods_curv) do method
        Ns, errs, errs_max = test_convergence(method, CurvatureScalar())
        method => scatterlines!(ax_curv, Ns, errs; label = string(method))
    end
    map(methods_tors) do method
        Ns, errs, errs_max = test_convergence(method, TorsionScalar())
        i = findfirst(splat((k, l) -> k === method), ls_curv)  # match colour of curvature line
        color = ls_curv[i][2].color
        scatterlines!(ax_tors, Ns, errs; label = string(method), color)
    end
    axislegend(ax_curv; position = (0.04, 0.35))
    let Ns = 2.0.^(range(7, 10; length = 3))
        let ys = @. (3 / Ns)^1
            lines!(ax_curv, Ns, ys; color = :grey, linestyle = :dash)
            text!(ax_curv, Ns[2], ys[2]; text = L"∼N^{-1}", color = :grey, fontsize = 20, align = (:left, :bottom))
        end
        let ys = @. (4 / Ns)^2
            lines!(ax_curv, Ns, ys; color = :grey, linestyle = :dash)
            text!(ax_curv, Ns[2], 0.8 * ys[2]; text = L"∼N^{-2}", color = :grey, fontsize = 20, align = (:right, :top))
        end
        let ys = @. (6 / Ns)^4
            lines!(ax_curv, Ns, ys; color = :grey, linestyle = :dash)
            text!(ax_curv, Ns[2], 0.8 * ys[2]; text = L"∼N^{-4}", color = :grey, fontsize = 20, align = (:right, :top))
        end
        let ys = @. (Ns / 1e9)^2
            lines!(ax_curv, Ns, ys; color = :grey, linestyle = :dash)
            text!(ax_curv, Ns[2], 0.8 * ys[2]; text = L"∼N^{2}", color = :grey, fontsize = 20, align = (:left, :top))
        end
        let ys = @. (8 / Ns)^4
            lines!(ax_tors, Ns, ys; color = :grey, linestyle = :dash)
            text!(ax_tors, Ns[2], 0.8 * ys[2]; text = L"∼N^{-4}", color = :grey, fontsize = 20, align = (:right, :top))
        end
        let ys = @. (Ns / 2e6)^3
            lines!(ax_tors, Ns, ys; color = :grey, linestyle = :dash)
            text!(ax_tors, Ns[2], 0.8 * ys[2]; text = L"∼N^{3}", color = :grey, fontsize = 20, align = (:left, :top))
        end
    end
    fig
end

if @isdefined(Makie)
    plot_convergence()
end
