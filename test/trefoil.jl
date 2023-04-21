using Test
using LinearAlgebra: norm, ⋅
using VortexPasta.Filaments
using VortexPasta.BiotSavart

function init_trefoil_filament(N::Int)
    R = π / 4
    S(t) = π .+ R .* Vec3(
        sinpi(t) + 2 * sinpi(2t),
        cospi(t) - 2 * cospi(2t),
        -sinpi(3t),
    )
    f = Filaments.init(ClosedFilament, N, CubicSplineMethod())
    ζs = range(0, 2; length = N + 1)[1:N]
    f .= S.(ζs)
    update_coefficients!(f)
    f
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

    cache_exact = @inferred(BiotSavart.init_cache(params_exact)).longrange
    cache_default = @inferred(BiotSavart.init_cache(params_default)).longrange

    @test BiotSavart.backend(cache_exact) isa ExactSumBackend
    @test BiotSavart.backend(cache_default) isa FINUFFTBackend

    # Compute induced velocity field in Fourier space
    BiotSavart.long_range_velocity_fourier!(cache_exact, fs)
    BiotSavart.long_range_velocity_fourier!(cache_default, fs)

    # Compare velocities in Fourier space.
    # Note: the comparison is not straightforward since the wavenumbers are not the same.
    # The "exact" implementation takes advantage of Hermitian symmetry to avoid
    # computing half of the data, while FINUFFT doesn't make this possible...
    ks_default = cache_default.wavenumbers
    ks_exact = cache_exact.wavenumbers
    inds_to_compare = ntuple(Val(3)) do i
        inds = eachindex(ks_exact[i])
        js = i == 1 ? inds[begin:end - 1] : inds
        @assert @views ks_default[i][js] == ks_exact[i][js]  # wavenumbers match in this index region
        js
    end |> CartesianIndices
    diffnorm_L2 = sum(inds_to_compare) do I
        uhat, vhat = cache_exact.uhat, cache_default.uhat
        du = uhat[I] - vhat[I]
        sum(abs2, du)  # = |u - v|^2
    end
    @test diffnorm_L2 < tol^2

    # Interpolate velocity back to filament positions
    BiotSavart.long_range_velocity_physical!(cache_exact, fs)
    BiotSavart.long_range_velocity_physical!(cache_default, fs)

    max_rel_error_physical = maximum(zip(cache_exact.charges, cache_default.charges)) do (qexact, qdefault)
        norm(qexact - qdefault) / norm(qexact)
    end
    @test max_rel_error_physical < tol

    # Copy data to arrays.
    vs_exact = map(f -> zero(Filaments.points(f)), fs)
    vs_default = map(f -> zero(Filaments.points(f)), fs)
    BiotSavart.add_long_range_velocity!(vs_exact, cache_exact)
    BiotSavart.add_long_range_velocity!(vs_default, cache_default)

    # Compare velocities one filament at a time.
    @test all(zip(vs_exact, vs_default)) do (u, v)
        isapprox(u, v; rtol = tol)
    end

    nothing
end

function compute_filament_velocity(f, α; params_kws...)
    params = ParamsBiotSavart(; params_kws..., α, rcut = 4 / α)
    cache = BiotSavart.init_cache(params)
    vs = zero(Filaments.points(f))
    BiotSavart.add_short_range_velocity_self!(vs, cache.shortrange, f)
    fs = [f]
    BiotSavart.long_range_velocity_fourier!(cache.longrange, fs)
    BiotSavart.add_long_range_velocity!(vs, cache.longrange, fs)
    vs
end

# Check that the total induced velocity doesn't depend strongly on the Ewald parameter α.
# (In theory it shouldn't depend at all...)
function check_independence_on_ewald_parameter(f, αs; params_kws...)
    vs_all = map(αs) do α
        compute_filament_velocity(
            f, α;
            # Use high-order quadratures to make sure that errors don't come from there.
            quadrature_short = GaussLegendreQuadrature(6),
            quadrature_long = GaussLegendreQuadrature(6),
            params_kws...,
        )
    end
    vs_test = first(vs_all)
    @test all(vs_all) do vs
        maxdiff = maximum(zip(vs, vs_test)) do (a, b)
            norm(a - b) / norm(b)
        end
        maxdiff < 1e-4
    end
    nothing
end

@testset "Trefoil" begin
    f = @inferred init_trefoil_filament(30)
    Ls = (2π, 2π, 2π)  # TODO test other sizes?
    Ns = (64, 64, 64)
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    params_kws = (; Ls, Ns, Γ = 2.0, a = 1e-5, α = kmax / 6,)
    @testset "Long range" begin
        compare_long_range([f]; tol = 1e-8, params_kws...)
    end
    @testset "Dependence on α" begin
        αs = [kmax / 5, kmax / 8, kmax / 16]
        check_independence_on_ewald_parameter(f, αs; params_kws...)
    end
end
