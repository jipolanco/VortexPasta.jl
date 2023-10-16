using Test
using LinearAlgebra: norm, ⋅, ×
using VortexPasta.PredefinedCurves: define_curve, TrefoilKnot
using VortexPasta.Filaments
using VortexPasta.BiotSavart

function init_trefoil_filament(N::Int; method = CubicSplineMethod())
    R = π / 4
    S = @inferred define_curve(TrefoilKnot(); translate = R, scale = R)
    ζs = range(0, 1; length = N + 1)[1:N]
    Xs = @inferred broadcast(S, ζs)  # same as S.(ζs)
    Filaments.init(ClosedFilament, Xs, method)
end

function compare_long_range(fs::AbstractVector{<:AbstractFilament}; tol = 1e-8, params_kws...)
    params_exact = @inferred ParamsBiotSavart(;
        params_kws...,
        backend_long = ExactSumBackend(),
    )
    params_default = @inferred ParamsBiotSavart(;
        params_kws...,
        backend_long = FINUFFTBackend(; tol,),
    )

    # This is just to check that Base.show is implemented for ParamsBiotSavart.
    @test startswith(repr(params_exact), "ParamsBiotSavart with:\n")
    @test startswith(repr(params_default), "ParamsBiotSavart with:\n")

    cache_exact_base = @inferred(BiotSavart.init_cache(params_exact, fs))
    cache_default_base = @inferred(BiotSavart.init_cache(params_default, fs))

    cache_exact = cache_exact_base.longrange
    cache_default = cache_default_base.longrange

    @test BiotSavart.backend(cache_exact) isa ExactSumBackend
    @test BiotSavart.backend(cache_default) isa FINUFFTBackend

    # Compute induced velocity field in Fourier space
    foreach((cache_exact, cache_default)) do c
        BiotSavart.compute_vorticity_fourier!(c, fs)
        BiotSavart.to_smoothed_velocity!(c)
    end

    # Compare velocities in Fourier space.
    # Note: the comparison is not straightforward since the wavenumbers are not the same.
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

    max_rel_error_physical = maximum(zip(cache_exact.common.charges, cache_default.common.charges)) do (qexact, qdefault)
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
        isapprox(u, v; rtol = tol)
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

    cache_naive = @inferred(BiotSavart.init_cache(params_naive, fs)).shortrange
    cache_cl = @inferred(BiotSavart.init_cache(params_cl, fs)).shortrange

    BiotSavart.set_filaments!(cache_naive, fs)
    BiotSavart.set_filaments!(cache_cl, fs)

    vs_naive = map(f -> zero(nodes(f)), fs)
    vs_cl = map(f -> zero(nodes(f)), fs)

    for (v, f) ∈ zip(vs_naive, fs)
        BiotSavart.add_short_range_velocity!(v, cache_naive, f)
    end

    for (v, f) ∈ zip(vs_cl, fs)
        BiotSavart.add_short_range_velocity!(v, cache_cl, f)
    end

    for (a, b) ∈ zip(vs_naive, vs_cl)
        @test isapprox(a, b; rtol = 1e-7)
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
            quadrature_short = quad,
            quadrature_long = quad,
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
        compare_long_range([f]; tol = 1e-8, params_kws..., α = kmax / 6)
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
end

##

if @isdefined(Makie)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    wireframe!(ax, Rect(0, 0, 0, Ls...); color = :grey, linewidth = 0.5)
    plot!(ax, f; refinement = 8)
    fig
end
