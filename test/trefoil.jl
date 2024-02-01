using Test
using LinearAlgebra: norm, ⋅, ×
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot
using VortexPasta.Filaments
using VortexPasta.Filaments: GeometricQuantity, Vec3
using VortexPasta.BiotSavart
using ForwardDiff: ForwardDiff
using LaTeXStrings  # used for plots only (L"...")

function trefoil_function()
    R = π / 4
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
            @test err < 1e-12
            @test err_max < 1e-11
        elseif quantity === CurvatureScalar()
            # Case of CurvatureScalar
            @test err < 1e-13
            @test err_max < 1e-13
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
        fs::AbstractVector{<:AbstractFilament}, backend = FINUFFTBackend();
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

    # This is just to check that Base.show is implemented for ParamsBiotSavart.
    @test startswith(repr(params_exact), "ParamsBiotSavart with:\n")
    @test startswith(repr(params_default), "ParamsBiotSavart with:\n")

    cache_exact_base = @inferred(BiotSavart.init_cache(params_exact, fs))
    cache_default_base = @inferred(BiotSavart.init_cache(params_default, fs))

    cache_exact = cache_exact_base.longrange
    cache_default = cache_default_base.longrange

    @test BiotSavart.backend(cache_exact) isa ExactSumBackend
    @test BiotSavart.backend(cache_default) === backend

    # Compute induced velocity field in Fourier space
    foreach((cache_exact, cache_default)) do c
        BiotSavart.add_point_charges!(c.common.pointdata, fs, quad)
        BiotSavart.compute_vorticity_fourier!(c)
        BiotSavart.to_smoothed_velocity!(c)
    end

    # Compare velocities in Fourier space.
    # Note: when using FINUFFTBackend, the comparison is not straightforward since the
    # wavenumbers are not the same.
    # The "exact" implementation takes advantage of Hermitian symmetry to avoid
    # computing half of the data, while FINUFFT doesn't make this possible...
    ks_default = cache_default.common.wavenumbers
    ks_exact = cache_exact.common.wavenumbers
    inds_to_compare = ntuple(Val(3)) do i
        inds = eachindex(ks_exact[i])
        js = i == 1 ? inds[begin:end - 1] : inds
        @assert @views ks_default[i][js] == ks_exact[i][js]  # wavenumbers match in this index region
        js
    end |> CartesianIndices
    diffnorm_L2 = sum(inds_to_compare) do I
        uhat, vhat = cache_exact.common.uhat, cache_default.common.uhat
        du = uhat[I] - vhat[I]
        sum(abs2, du)  # = |u - v|^2
    end
    @test diffnorm_L2 < tol^2

    # Interpolate velocity back to filament positions
    foreach((cache_exact, cache_default)) do c
        BiotSavart.set_interpolation_points!(c, fs)
        BiotSavart.interpolate_to_physical!(c)
    end

    charges(cache) = cache.common.pointdata.charges  # get charges from cache
    max_rel_error_physical = maximum(zip(charges(cache_exact), charges(cache_default))) do (qexact, qdefault)
        norm(qexact - qdefault) / norm(qexact)
    end
    @test max_rel_error_physical < tol

    # Copy data to arrays.
    vs_exact = map(f -> zero(nodes(f)), fs)
    vs_default = map(f -> zero(nodes(f)), fs)
    BiotSavart.add_long_range_output!(vs_exact, cache_exact)
    BiotSavart.add_long_range_output!(vs_default, cache_default)

    # Compare velocities one filament at a time.
    check_approx(us, vs, tol) = all(zip(us, vs)) do (u, v)
        # Note: we avoid comparing ghost cells in case these are PaddedVectors.
        @views isapprox(u[:], v[:]; rtol = tol)
    end
    @test check_approx(vs_default, vs_exact, tol)

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

    @test check_approx(vs_default, vs_exact, tol)
    @test check_approx(ψs_default, ψs_exact, tol)

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

    cache_naive = @inferred BiotSavart.init_cache(params_naive, fs)
    cache_cl = @inferred BiotSavart.init_cache(params_cl, fs)

    vs_naive = map(similar ∘ nodes, fs)
    vs_cl = map(similar ∘ nodes, fs)

    BiotSavart.velocity_on_nodes!(vs_naive, cache_naive, fs; longrange = Val(false))
    BiotSavart.velocity_on_nodes!(vs_cl, cache_cl, fs; longrange = Val(false))

    for (a, b) ∈ zip(vs_naive, vs_cl)
        # Note: since `a` and `b` are PaddedVectors, we need to exclude their ghost cells
        # (which can be different) from the comparison. For that reason we use a view onto
        # the "central" (non-ghost) data of both arrays.
        @test @views isapprox(a[:], b[:]; rtol = 1e-7)
    end

    nothing
end

function compute_filament_velocity_and_streamfunction(f::AbstractFilament; α, Ls, params_kws...)
    rcut = 4 * sqrt(2) / α
    @assert rcut < minimum(Ls) / 3  # cell lists requirement
    params = ParamsBiotSavart(; params_kws..., α, Ls, rcut)
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
function check_independence_on_ewald_parameter(f, αs; quad = GaussLegendre(2), params_kws...)
    fields_all = map(αs) do α
        compute_filament_velocity_and_streamfunction(
            f;
            α,
            backend_short = CellListsBackend(2),
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
    @show maxdiffs_vel maxdiffs_stf
    @test maximum(maxdiffs_vel) < 1e-4
    @test maximum(maxdiffs_stf) < 1e-5
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
    f = @inferred init_trefoil_filament(30)
    Ls = (1.5π, 1.5π, 2π)  # Ly is small to test periodicity effects
    Ns = (3, 3, 4) .* 30
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    params_kws = (; Ls, Ns, Γ = 2.0, a = 1e-5,)
    @testset "Long range" begin
        # Test NUFFT backends with default parameters
        compare_long_range([f], NonuniformFFTsBackend(); tol = 1e-8, params_kws..., α = kmax / 6)
        compare_long_range([f], FINUFFTBackend(); tol = 1e-8, params_kws..., α = kmax / 6)
        @testset "$quantity near r = 0" for quantity ∈ (Velocity(), Streamfunction())
            test_long_range_accuracy_near_zero(Float64, quantity)
        end
    end
    @testset "Short range" begin
        compare_short_range([f]; params_kws..., α = kmax / 6)
    end
    @testset "Dependence on α" begin
        αs = [kmax / 5, kmax / 6, kmax / 7, kmax / 8, kmax / 12, kmax / 16]
        quadratures = (GaussLegendre(3), NoQuadrature())
        @testset "$quad" for quad ∈ quadratures
            check_independence_on_ewald_parameter(f, αs; quad, params_kws...)
        end
        @testset "FourierMethod()" begin
            f_fourier = @inferred init_trefoil_filament(32; method = FourierMethod())
            check_independence_on_ewald_parameter(f_fourier, αs; params_kws...)
        end
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
