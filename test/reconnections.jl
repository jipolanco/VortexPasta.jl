using Test
using VortexPasta.PredefinedCurves:
    define_curve, Ring, TrefoilKnot, Lemniscate, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3
using VortexPasta.FilamentIO
using VortexPasta.Reconnections
using VortexPasta.Reconnections: reconnect!
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Timestepping: VortexFilamentSolver
using VortexPasta.Diagnostics
using StaticArrays
using Rotations
using UnicodePlots: lineplot, lineplot!
using LinearAlgebra: norm, I, ⋅
using JET: JET
using KernelAbstractions: KernelAbstractions as KA  # for JET only

VERBOSE::Bool = get(ENV, "JULIA_TESTS_VERBOSE", "false") in ("true", "1")

# This is just to fake the `info` argument of `reconnect_with_itself!` and
# `reconnect_with_other!`, which is usually returned by the `should_reconnect` function.
function dummy_reconnection_info(::Type{T}) where {T}
    (;
        p⃗ = zero(Vec3{T}),
        length_before = zero(T),
        length_after = zero(T),
    )
end

# Create a curve which resembles an 8 (or ∞).
# This specific curve is called the lemniscate of Bernoulli (https://en.wikipedia.org/wiki/Lemniscate_of_Bernoulli)
# We allow a perturbation in the 3rd direction (Az) so that the curve doesn't exactly cross itself.
function figure_eight_curve(; a::Real = 1, translate = π * Vec3(1, 1, 1), Az = 0)
    # Note: the crossings are located at t = 0.25 and 0.75.
    tlims = (0, 1)
    S = define_curve(Lemniscate(; Az,); translate, scale = a,)
    (; tlims, S,)
end

function ring_curve(; scale = 1, translate = π * Vec3(1, 1, 1), sign = +1, rotate = I)
    tlims = (0, 1)
    S = define_curve(Ring(); translate, scale, orientation = sign, rotate)
    (; tlims, S,)
end

# Here `orientation ∈ 1:3` refers to the direction of the line (x, y or z).
# Not to be confused with the `orientation` argument of `define_curve`.
function infinite_line_curve(; origin, L = 2π, sign, orientation::Int)
    tlims = (0, 1)
    rotate = if orientation == 1
        RotMatrix(RotY(π / 2))
    elseif orientation == 2
        RotMatrix(RotX(-π / 2))
    elseif orientation == 3
        RotMatrix(RotZ(0))  # identity matrix (no rotation)
    end
    S = define_curve(
        PeriodicLine();
        scale = L, orientation = sign, rotate, translate = origin,
    )
    offset = S(1) - S(0)
    offset_alt = setindex(zero(Vec3{Float64}), sign * L, orientation)
    @assert offset == offset_alt
    (; tlims, S, offset,)
end

# Check that the filament doesn't have any crazy jumps between two nodes.
function no_jumps(f::AbstractFilament, lmax)
    Xs = parent(nodes(f))  # use `parent` to consider periodically-padded nodes as well
    maximum(norm, diff(Xs)) < lmax
end

function compare_with_filament(f::AbstractFilament, periods; kwargs...)
    to_main_cell(x⃗) = mod.(x⃗, periods)
    function (g::AbstractFilament)
        length(f) == length(g) || return false
        # Try to find offset between the two filaments.
        x = to_main_cell(f[begin])
        noff = firstindex(g)
        while noff <= lastindex(g) && to_main_cell(g[noff]) != x
            noff += 1
        end
        noff > lastindex(g) && return false
        noff -= 1
        for i in eachindex(f)
            j = mod1(i + noff, length(g))
            x = to_main_cell(f[i])
            y = to_main_cell(g[j])
            isapprox(x, y; kwargs...) || return false
        end
        true
    end
end

function write_vectors(io, vecs...)
    for i in eachindex(vecs...)
        vals = map(v -> v[i], vecs)
        join(io, vals, '\t')
        print(io, '\n')
    end
end

trefoil_scheme_dt(::Any) = trefoil_scheme_dt(RK4())  # default
trefoil_scheme_dt(::RK4) = 1.0
trefoil_scheme_dt(::DP5) = 1.0
trefoil_scheme_dt(::SSPRK33) = 0.8
trefoil_scheme_dt(::KenCarp3) = 0.5
trefoil_scheme_dt(::Midpoint) = 0.3

# The `large_trefoil = true` case was added after the "restart" bug (not really related to
# restarts) was fixed. It allows to check that filament buffers stay valid after
# reconnections involving a change of end-to-end offset.
function test_trefoil_knot_reconnection(
        ::Type{Criterion} = ReconnectBasedOnDistance;
        large_trefoil = false,
        scheme = RK4(),
    ) where {Criterion <: ReconnectionCriterion}
    test_jet = get(ENV, "JULIA_ENABLE_JET_KA_TESTS", "false") ∈ ("true", "1")  # enable JET tests involving KA kernels
    scale = large_trefoil ? (π / 2.8) : (π / 4)
    S = define_curve(TrefoilKnot(); translate = π, scale)
    N = 64
    f = Filaments.init(S, ClosedFilament, N, QuinticSplineMethod())
    fs_init = [f]

    params_bs = let L = 2π
        Ls = (1, 1, 1) .* L
        β = 3.5  # accuracy coefficient
        rcut = L / 4
        α = β / rcut
        kmax = 2α * β
        M = ceil(Int, kmax * L / π)
        Ns = (1, 1, 1) .* M
        ParamsBiotSavart(;
            Γ = 2.0, α, a = 1e-6, Δ = 1/4, rcut, Ls, Ns,
            backend_short = CellListsBackend(2),
            backend_long = NonuniformFFTsBackend(σ = 1.5, m = HalfSupport(4)),
            quadrature = GaussLegendre(3),
            quadrature_near_singularity = GaussLegendre(3),
            lia_segment_fraction = 0.2,
        )
    end
    # println(params_bs)

    times = Float64[]
    energy = similar(times)
    helicity = similar(times)
    torsion_integral = similar(times)
    n_reconnect = Ref(0)
    t_reconnect = Ref(0.0)
    n_reconnections = Ref(0)

    function callback(iter)
        local (; fs, t, dt, nstep,) = iter
        local (; Γ,) = iter.prob.p
        if nstep == 0
            foreach(empty!, (times, energy))
            n_reconnect[] = t_reconnect[] = 0
        end
        push!(times, t)
        local quad = GaussLegendre(4)
        local E = Diagnostics.kinetic_energy(iter; quad)
        local H = Diagnostics.helicity(iter; quad)
        push!(energy, E)
        push!(helicity, H)
        Nf = length(fs)

        # Compute torsion integral
        local tor = zero(eltype(torsion_integral))
        for f in fs
            tor += integrate(f, quad) do ff, i, ζ
                local s⃗′ = ff(i, ζ, Derivative(1))
                local τ = ff(i, ζ, TorsionScalar())
                τ * norm(s⃗′)
            end
        end
        tor /= 2π
        push!(torsion_integral, tor)

        if iter.stats.reconnection_count > n_reconnections[]
            rec_before = n_reconnections[]
            rec_after = iter.stats.reconnection_count
            n_reconnections[] = rec_after
            if VERBOSE
                @info "Performed $(rec_after - rec_before) reconnections" rec_before rec_after
            end
        end

        # save_checkpoint("trefoil_$nstep.vtkhdf", iter; refinement = 4, periods = nothing) do io
        #     io["velocity"] = iter.vs
        #     io["streamfunction"] = iter.ψs
        # end
        # save_checkpoint("trefoil_points_$nstep.vtkhdf", iter; periods = nothing, refinement = 1) do io
        #     io["velocity"] = iter.vs
        #     io["streamfunction"] = iter.ψs
        # end
        if VERBOSE
            @show nstep, t, Nf, E, H/Γ^2, tor
        end
        if n_reconnect[] == 0 && Nf == 2
            t_reconnect[] = t
            n_reconnect[] = nstep
        end
        nothing
    end

    tspan = (0.0, 2.0)
    prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs)
    δ = Filaments.minimum_node_distance(prob.fs)
    if Criterion <: ReconnectBasedOnDistance
        d_crit = 0.75 * δ
        reconnect = ReconnectBasedOnDistance(d_crit; max_passes = 1, use_velocity = true)
    elseif Criterion <: ReconnectFast
        d_crit = 0.75 * δ
        reconnect = ReconnectFast(d_crit; nthreads = 1, max_passes = 4, use_velocity = true)  # we use nthreads = 1 to test the serial implementation
    end
    dt = BiotSavart.kelvin_wave_period(params_bs, δ) * trefoil_scheme_dt(scheme)
    # @show δ d_crit dt
    iter = @inferred init(
        prob, scheme;
        dt,
        refinement = RefineBasedOnSegmentLength(0.75 * δ),
        # refinement = RefineBasedOnCurvature(0.4; ℓ_max = 1.5 * l_min, ℓ_min = 0.4 * l_min),
        reconnect,
        adaptivity = NoAdaptivity(),
        filament_nderivs = Val(3),  # for computation of torsion
        callback,
    )

    let
        ks, Ek = Diagnostics.energy_spectrum(iter)
        _, Hk = Diagnostics.helicity_spectrum(iter)
        # open(io -> write_vectors(io, ks, Ek, Hk), "trefoil_spectra_before.dat", "w")
    end

    if VERBOSE
        @time solve!(iter)
        println(iter.to)
    else
        solve!(iter)
    end

    for f in iter.fs
        offset = Filaments.end_to_end_offset(f)
        if large_trefoil
            # End-to-end offsets are nonzero in "large" case (this is the whole purpose of
            # the test).
            @test offset != zero(offset)
        else
            @test offset == zero(offset)  # end-to-end offsets are zero in "small" case
        end
    end

    let
        ks, Ek = Diagnostics.energy_spectrum(iter)
        _, Hk = Diagnostics.helicity_spectrum(iter)
        # open(io -> write_vectors(io, ks, Ek, Hk), "trefoil_spectra_after.dat", "w")
    end

    energy_rel = energy ./ energy[begin]

    if VERBOSE
        plt = lineplot(
            times, energy_rel;
            xlabel = "Time", ylabel = "Energy", title = "Trefoil knot / $scheme",
        )
        println(plt)
        plt = lineplot(
            times, helicity;
            xlabel = "Time", ylabel = "Helicity", title = "Trefoil knot / $scheme",
        )
        # lineplot!(plt, times, torsion_integral)
        println(plt)
    end

    if VERBOSE
        @show t_reconnect[]
        @show last(energy_rel)
    end

    let stats = iter.stats
        VERBOSE && @show stats
        @test stats.reconnection_count == 3        # total number of reconnections
        @test stats.reconnection_length_loss > 0   # loss of vortex length due to reconnections
        @test stats.filaments_removed_count == 0   # no filaments were removed
        @test stats.filaments_removed_length == 0  # no filaments were removed
    end

    @test n_reconnect[] > 0
    if large_trefoil
        @test 0.57 < t_reconnect[] < 0.67
        @test 0.965 < last(energy_rel) < 0.985
    elseif Criterion <: ReconnectBasedOnDistance
        @test 1.60 < t_reconnect[] < 1.70  # this depends on several parameters...
        @test 0.98 < last(energy_rel) < 0.99
    elseif Criterion <: ReconnectFast
        @test 1.60 < t_reconnect[] < 1.70
        @test 0.970 < last(energy_rel) < 0.980
    end

    if test_jet
        JET.@test_opt VortexFilamentProblem(fs_init, tspan, params_bs)
        JET.@test_opt ignored_modules=(Base, KA, Base.IteratorsMD) step!(iter)
    end

    nothing
end

function test_static_figure_eight_knot(
        ::Type{Criterion} = ReconnectBasedOnDistance;
    ) where {Criterion <: ReconnectionCriterion}
    Az = 0.001
    curve = figure_eight_curve(; a = π / 4, Az, translate = (0, 0, 0))
    (; S, tlims,) = curve
    N = 32
    ζs = range(tlims...; length = 2N + 1)[2:2:2N]
    f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())

    l_min = minimum_knot_increment(f)

    if Criterion <: ReconnectBasedOnDistance
        # This doesn't really depend on the reconnection criterion
        @testset "Filaments.split!" begin
            i = length(f) ÷ 4
            j = 3i
            JET.@test_opt ignored_modules=(Base,) Filaments.split!(copy(f), i, j)
            f1, f2 = @inferred Filaments.split!(copy(f), i, j)
            update_coefficients!.((f1, f2))
            @test length(f1) == length(f2) == length(f) ÷ 2
            @test f1 == f[i + 1:j]
            @test f2 == vcat(f[j + 1:end], f[begin:i])
        end
    end

    @testset "reconnect!" begin
        if Criterion <: ReconnectBasedOnDistance
            crit = @inferred ReconnectBasedOnDistance(l_min / 2; max_passes = 2, use_velocity = false)
        elseif Criterion <: ReconnectFast
            crit = @inferred ReconnectFast(0.9 * l_min; use_velocity = false)
        end
        fc = copy(f)
        fs_all = [fc]
        Ls = (2π, 2π, 2π)
        cache = @inferred Reconnections.init_cache(crit, fs_all, Ls)
        # @test_opt ignored_modules=(Base,) reconnect!(cache, fs_all)
        rec = @inferred reconnect!(cache, fs_all)
        @test rec.reconnection_count == 1
        @test length(fs_all) == 2  # filament split into two!
        f1, f2 = fs_all[1], fs_all[2]

        # # Check that the reconnection happened at the right place
        # i = length(f) ÷ 4
        # j = 3i
        # @test length(f1) == length(f2) == length(f) ÷ 2
        # @test f1 == f[i + 1:j]
        # @test f2 == vcat(f[j + 1:end], f[begin:i])
    end
end

function test_static_antiparallel_rings(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    # Test two nearly overlapping rings with the same sign (they should reconnect,
    # since they're antiparallel at the reconnection point).
    R = π / 4
    ϵ = 0.02R
    rings = (
        ring_curve(; scale = R, translate = Vec3(π - (R + ϵ), π, π), sign = +1),
        ring_curve(; scale = R, translate = Vec3(π + (R + ϵ), π, π), sign = +1),
    )
    d_min = norm(rings[2].S(0.5) - rings[1].S(0.0))  # minimum distance between the rings
    @assert d_min ≈ 2ϵ

    fs_orig = map(rings) do ring
        (; S, tlims,) = ring
        N = 32
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    end
    l_min = minimum(minimum_knot_increment, fs_orig)
    d_min_nodes = let (f, g) = fs_orig
        minimum(Iterators.product(eachindex(f), eachindex(g))) do (i, j)
            norm(f[i] - g[j])
        end
    end
    if Criterion <: ReconnectBasedOnDistance
        crit = ReconnectBasedOnDistance(1.2 * d_min_nodes; use_velocity = false)
    elseif Criterion <: ReconnectFast
        crit = ReconnectFast(1.2 * d_min_nodes; use_velocity = false)
    end

    @testset "reconnect!" begin
        fs = collect(copy.(fs_orig))
        if Criterion <: ReconnectFast
            Ls = (2π, 2π, 2π)  # ReconnectFast requires periodic domain
            cache = @inferred Reconnections.init_cache(crit, fs, Ls)
        else
            cache = @inferred Reconnections.init_cache(crit, fs)
        end
        rec = @inferred reconnect!(cache, fs)
        @test rec.reconnection_count == 1  # one reconnection
        @test length(fs) == 1
        @test length(fs[1]) == sum(length, fs_orig)
        @test maximum_knot_increment(fs) < 2 * l_min  # check that there are no crazy jumps
    end

    nothing
end

function test_static_nearly_overlapping_rings(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    # Test two nearly overlapping rings (with some small offset in z) which should reconnect at two points.
    # We should end up with two separate vortices.
    # With ReconnectBasedOnDistance, this requires two reconnection steps:
    #   (1) the two vortices merge at one of the reconnection points, and
    #   (2) the resulting vortex splits into two at the other reconnection point.
    R = π / 4
    ϵ = 0.02R
    rings = (
        ring_curve(; scale = R, translate = Vec3(π - R / sqrt(2), π, π + ϵ), sign = +1),
        ring_curve(; scale = R, translate = Vec3(π + R / sqrt(2), π, π - ϵ), sign = +1),
    )

    fs_orig = map(rings) do ring
        (; S, tlims,) = ring
        N = 32
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    end
    d_min_nodes = let (f, g) = fs_orig
        minimum(Iterators.product(eachindex(f), eachindex(g))) do (i, j)
            norm(f[i] - g[j])
        end
    end
    fs = collect(copy.(fs_orig))

    if Criterion <: ReconnectBasedOnDistance
        crit = ReconnectBasedOnDistance(1.2 * d_min_nodes; use_velocity = false)
        cache = @inferred Reconnections.init_cache(crit, fs)
        # We need two reconnect! passes to arrive to the final state.
        for n ∈ 1:4
            rec = reconnect!(cache, fs)
            nrec = rec.reconnection_count
            if n ≤ 2
                @test nrec == 1
            else
                @test nrec == 0
                break
            end
        end
    elseif Criterion <: ReconnectFast
        crit = ReconnectFast(1.2 * d_min_nodes; use_velocity = false)
        Ls = (2π, 2π, 2π)  # ReconnectFast requires periodic domain
        cache = @inferred Reconnections.init_cache(crit, fs, Ls)
        # This criterion should perform all reconnections in a single pass.
        rec = reconnect!(cache, fs)
        @test rec.reconnection_count == 2
    end

    @test length(fs) == 2
    @test fs[1] != fs_orig[1]  # reconnection happened
    @test fs[2] != fs_orig[2]  # reconnection happened
    @test length(fs[1]) ≠ length(fs[2])  # the resulting vortices have different sizes

    nothing
end

function test_static_periodic_nearly_overlapping_ellipses(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    # Periodic ellipses which almost touch (distance 2ϵ).
    # A single ellipse reconnects into 2 infinite lines.
    # We use an ellipse instead of a circle to make sure that reconnections happen
    # along one direction only (simpler).
    # Note that we could also use a circular ring with different domain periodicities...
    periods = 2π .* (1, 1, 1)
    R = π
    ϵ = π / 100
    N = 32
    scale = R * SDiagonal(1, 0.5, 1)  # "squeeze" ring along `y` direction
    ring = ring_curve(; scale, translate = π * Vec3(1, 1, 1),)
    f_orig = let
        (; S, tlims,) = ring
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    end

    l_min = minimum_knot_increment(f_orig)

    fs_orig = [copy(f_orig)]
    crit = Criterion(4 * ϵ; use_velocity = false)
    fs = copy.(fs_orig)
    cache = @inferred Reconnections.init_cache(crit, fs, periods)
    rec = @inferred reconnect!(cache, fs)

    @test rec.reconnection_count == 1
    @test length(fs) == 2
    @test abs(end_to_end_offset(fs[1])[1]) ≈ 2π  # infinite line in x
    @test abs(end_to_end_offset(fs[2])[1]) ≈ 2π  # infinite line in x
    foreach(Filaments.update_coefficients!, fs)
    @test all(f -> no_jumps(f, 2 * l_min), fs)

    nothing
end

function test_static_periodic_overlapping_ellipses(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    # A single ellipse whose major axis is slightly larger than the domain period.
    # The ellipse reconnects onto 2 infinite lines (similar to previous test) plus a
    # small ring.
    periods = 2π .* (1, 1, 1)
    R = π * 1.2
    N = 64
    scale = R * SDiagonal(1, 0.5, 1)  # "squeeze" ring along `y` direction
    ring = ring_curve(; scale, translate = π * Vec3(1, 1, 1),)
    f_orig = let
        (; S, tlims,) = ring
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), FiniteDiffMethod(2))
    end

    l_min = minimum_knot_increment(f_orig)
    fs = [copy(f_orig)]
    if Criterion <: ReconnectBasedOnDistance
        crit = ReconnectBasedOnDistance(0.5 * l_min; use_velocity = false)
    elseif Criterion <: ReconnectFast
        crit = ReconnectFast(0.9 * l_min; use_velocity = false)
    end
    cache = @inferred Reconnections.init_cache(crit, fs, periods)

    if Criterion <: ReconnectBasedOnDistance
        # We need 2 passes to arrive to the final state.
        for n ∈ 1:4
            rec = @inferred reconnect!(cache, fs)
            nrec = rec.reconnection_count
            if n ≤ 2
                @test nrec == 1
            else
                @test nrec == 0
                break
            end
        end
    elseif Criterion <: ReconnectFast
        # We just need a single pass.
        rec = @inferred reconnect!(cache, fs)
        @test rec.reconnection_count == 2
    end

    @test length(fs) == 3
    @test sum(length, fs) == N  # number of nodes didn't change

    @test sum(f -> norm(end_to_end_offset(f)), fs) ≈ 4π  # two infinite lines
    @test all(f -> no_jumps(f, 2.1 * l_min), fs)
    # @show end_to_end_offset.(fs)

    # Check that we obtain 3 filaments with orientations -2π, +2π and 0 along x.
    offsets = end_to_end_offset.(fs)
    @test count(isapprox(Vec3(-2π, 0, 0)), offsets) == 1  # infinite line
    @test count(isapprox(Vec3(+2π, 0, 0)), offsets) == 1  # infinite line
    @test count(iszero, offsets) == 1  # small loop

    nothing
end

function test_static_periodic_overlapping_circles(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    # Slightly more complex case of overlapping periodic circles.
    # Each circle is reconnected into 3 closed curves.
    # With ReconnectBasedOnDistance, this requires 4 passes of reconnect!.
    Ls = 2π .* (1, 1, 1)
    R = π * 1.2
    N = 64
    ring = ring_curve(; scale = R, translate = π * Vec3(1, 1, 1),)
    f_orig = let
        (; S, tlims,) = ring
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), FiniteDiffMethod(2))
    end
    fs_orig = [f_orig]

    l_min = minimum_knot_increment(f_orig)
    if Criterion <: ReconnectBasedOnDistance
        crit = ReconnectBasedOnDistance(0.5 * l_min; use_velocity = false)
    elseif Criterion <: ReconnectFast
        crit = ReconnectFast(0.6 * l_min; use_velocity = false)
    end

    fs = copy.(fs_orig)
    cache = @inferred Reconnections.init_cache(crit, fs, Ls)

    if Criterion <: ReconnectBasedOnDistance
        nfilaments = (2, 3, 2, 3)  # expected number of filaments after each pass
        for n ∈ 1:10
            # write_vtkhdf("reconnect_big_circle_$(Criterion)_step$(n - 1).vtkhdf", fs)
            rec = @inferred reconnect!(cache, fs)
            nrec = rec.reconnection_count
            if n ≤ 4
                @test length(fs) == nfilaments[n]
                @test nrec == 1
            else
                @test nrec == 0
                break
            end
        end
    else
        rec = @inferred reconnect!(cache, fs)
        # write_vtkhdf("reconnect_big_circle_$(Criterion).vtkhdf", fs)
        @test rec.reconnection_count == 4
    end

    @test length(fs) == 3
    @test sum(length, fs) == N  # number of nodes didn't change
    @test all(f -> iszero(end_to_end_offset(f)), fs)  # all closed filaments
    @test all(f -> no_jumps(f, 2.1 * l_min), fs)

    if Criterion <: ReconnectBasedOnDistance
        @testset "Using max_passes = 10" begin
            crit = @inferred ReconnectBasedOnDistance(l_min / 2; max_passes = 10)
            fs = copy.(fs_orig)
            cache = @inferred Reconnections.init_cache(crit, fs, Ls)
            rec = @inferred reconnect!(cache, fs)
            @test rec.reconnection_count == 4  # total number of reconnections
            @test rec.npasses == 5  # total number of passes required (the last one doesn't do any reconnections)
            @test length(fs) == 3
            @test sum(length, fs) == N  # number of nodes didn't change
            @test all(f -> iszero(end_to_end_offset(f)), fs)  # all closed filaments
            @test all(f -> no_jumps(f, 2.1 * l_min), fs)
        end
    end

    nothing
end

function test_static_periodic_line_with_ring(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    R = 0.6π
    lines = [
        infinite_line_curve(; L = 2π, orientation = 1, origin = Vec3(1, 1, 1) * π, sign = +1),
        ring_curve(; scale = R, translate = Vec3(1, 1, 1) * π, sign = +1,)
    ]
    Ns = 10, 12
    fs_orig = map(lines, Ns) do line, N
        (; S, tlims,) = line
        offset = S(tlims[2]) - S(tlims[1])
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod(); offset)
    end
    fs = copy.(fs_orig)
    l_min = minimum_knot_increment(fs_orig)
    periods = 2π .* (1, 1, 1)

    # Manual reconnection: merge then split
    # 1. Determine merge location
    x_merge = π - R
    i = searchsortedlast(fs[1], x_merge; by = X -> X[1])
    _, j = findmin(X -> X[1], fs[2])
    @assert fs[1][i].x ≤ x_merge ≤ fs[1][i + 1].x
    let a = fs[2][j], b = fs[2][j + 1]
        @assert π - a.y ≈ b.y - π
        @assert a.x ≈ b.x
        @assert fs[1][i].x < a.x < fs[1][i + 1].x
    end

    h = Filaments.merge!(fs[1], fs[2], i, j)
    Filaments.update_coefficients!(h)
    @test end_to_end_offset(h) == end_to_end_offset(fs_orig[1])
    @test no_jumps(h, 2 * l_min)
    fs = [h]

    # 2. Determine split location
    x_split = π + R
    j = findlast(X -> X[1] < x_split, h)
    @assert h[j].x ≤ x_split < h[j + 1].x
    i = findlast(X -> X[1] > h[j].x, @view(h[begin:j - 1])) - 1
    let a = h[i], b = h[i + 1]
        @assert π - a.y ≈ b.y - π
        @assert a.x ≈ b.x
        @assert h[j].x < a.x < h[j + 1].x
    end

    f, g = Filaments.split!(h, i, j)
    fs_manual = [f, g]
    for f ∈ fs_manual
        Filaments.fold_periodic!(f, periods)
        Filaments.update_coefficients!(f)
    end
    @test iszero(end_to_end_offset(f))  # closed loop
    @test end_to_end_offset(g) == end_to_end_offset(fs_orig[1])  # infinite line
    @test no_jumps(f, 2 * l_min)
    @test no_jumps(g, 2 * l_min)

    # Automatic reconnection
    fs = copy.(fs_orig)
    if Criterion <: ReconnectBasedOnDistance
        crit = ReconnectBasedOnDistance(0.5 * l_min; use_velocity = false)
    elseif Criterion <: ReconnectFast
        crit = ReconnectFast(0.9 * l_min; use_velocity = false)
    end

    cache = @inferred Reconnections.init_cache(crit, fs, periods)

    if Criterion <: ReconnectBasedOnDistance
        for n ∈ 1:5
            rec = reconnect!(cache, fs)
            nrec = rec.reconnection_count
            if n ≤ 2
                @test nrec == 1
            else
                @test nrec == 0
                break
            end
        end
    else
        rec = reconnect!(cache, fs)
        @test rec.reconnection_count == 2
    end

    for f ∈ fs
        Filaments.fold_periodic!(f, periods)
        Filaments.update_coefficients!(f)
    end

    @test length(fs) == 2

    # Compare with manually merged vortices
    rtol = 1e-15
    @test count(compare_with_filament(fs_manual[1], periods; rtol), fs) == 1
    @test count(compare_with_filament(fs_manual[2], periods; rtol), fs) == 1

    nothing
end

function test_static_periodic_line_with_line(
        ::Type{Criterion};
    ) where {Criterion <: ReconnectionCriterion}
    lines = [
        infinite_line_curve(; L = 2π, orientation = 1, origin = Vec3(1, 1, 0.98) * π, sign = +1),
        infinite_line_curve(; L = 2π, orientation = 2, origin = Vec3(1, 1, 1.02) * π, sign = +1),
    ]
    Ns = 8, 12
    fs_orig = map(lines, Ns) do line, N
        (; S, offset, tlims,) = line
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod(); offset)
    end

    l_min = minimum_knot_increment(fs_orig)

    # Manual merge
    Xmerge = Vec3(1, 1, 1) * π
    fs = copy.(fs_orig)
    i = findlast(X -> X[1] < Xmerge[1], fs[1])
    j = findlast(X -> X[2] < Xmerge[2], fs[2])
    h = Filaments.merge!(fs[1], fs[2], i, j)
    Filaments.update_coefficients!(h)
    fs_manual = [h]
    @test no_jumps(h, 2 * l_min)
    @test end_to_end_offset(h) == sum(end_to_end_offset, fs_orig)  # output offset = sum of input offsets

    # Automatic reconnection
    crit = Criterion(l_min; use_velocity = false)
    periods = 2π .* (1, 1, 1)
    fs = copy.(fs_orig)
    cache = @inferred Reconnections.init_cache(crit, fs, periods)
    rec = reconnect!(cache, fs)
    @test rec.reconnection_count == 1
    @test length(fs) == 1

    # Compare with manually merged vortices
    @test compare_with_filament(fs_manual[1], periods; rtol = 1e-15)(fs[1])

    nothing
end

##

@testset "Reconnections" begin
    # Test filament removal after a self-reconnection, when one or the two resulting
    # filaments is too small (≤ 3 nodes for CubicSplineMethod).
    # Here we test the internal function `reconnect_with_itself!`.
    # (This is used in ReconnectBasedOnDistance only.)
    @testset "Static: remove filaments" begin
        Az = 0.001
        curve = figure_eight_curve(; a = π / 4, Az, translate = (0, 0, 0))
        (; S, tlims,) = curve
        @testset "Remove zero or one filament" begin
            N = 8  # very few points to ensure that a single filament will be removed
            ζs = range(tlims...; length = 2N + 1)[2:2:2N]
            info = dummy_reconnection_info(Float64)
            let i = 1, j = 5  # keep both filaments (both have 4 nodes)
                f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
                fs = [f]
                nremoved, Lremoved = @inferred Reconnections.reconnect_with_itself!(Returns(nothing), fs, f, i, j, info)
                @test nremoved == 0
                @test Lremoved == 0
                @test length(fs) == 2  # the new `fs` vector contains two filaments
            end
            let i = 1, j = 3  # remove the first filament (2 nodes)
                f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
                fs = [f]
                nremoved, Lremoved = @inferred Reconnections.reconnect_with_itself!(Returns(nothing), fs, f, i, j, info)
                @test nremoved == 1
                # @show Lremoved
                @test 0.88 < Lremoved < 0.89
                @test length(fs) == 1  # the new `fs` vector contains a single filament
            end
            let i = 1, j = 7  # remove the second filament (2 nodes)
                f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
                fs = [f]
                nremoved, Lremoved = @inferred Reconnections.reconnect_with_itself!(Returns(nothing), fs, f, i, j, info)
                @test nremoved == 1
                # @show Lremoved
                @test 0.96 < Lremoved < 0.97
                @test length(fs) == 1  # the new `fs` vector contains a single filament
            end
        end
        @testset "Remove the two filaments" begin
            N = 4  # very few points to ensure that both filaments will be removed
            ζs = range(tlims...; length = 2N + 1)[2:2:2N]
            info = dummy_reconnection_info(Float64)
            let i = 1, j = 3  # remove the two filaments (2 nodes each)
                f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
                fs = [f]
                nremoved, Lremoved = @inferred Reconnections.reconnect_with_itself!(Returns(nothing), fs, f, i, j, info)
                @test nremoved == 2
                # @show Lremoved
                @test 2.0 < Lremoved < 2.1
                @test length(fs) == 0  # the new `fs` vector is empty
            end
        end
    end

    @testset "Dynamic: trefoil knot" begin
        schemes = (RK4(),)
        criteria = (ReconnectBasedOnDistance, ReconnectFast)
        @testset "Scheme: $scheme" for scheme in schemes
            @testset "Criterion: $Criterion" for Criterion in criteria
                @testset "Large trefoil: $large_trefoil" for large_trefoil in (false, true)
                    test_trefoil_knot_reconnection(Criterion; scheme, large_trefoil)
                end
            end
        end
    end

    @testset "Static self-reconnection: figure-eight knot" begin
        for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
            @testset "$Criterion" test_static_figure_eight_knot(Criterion)
        end
    end

    @testset "Static: antiparallel rings" begin
        for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
            @testset "$Criterion" test_static_antiparallel_rings(Criterion)
        end
    end

    @testset "Static: nearly overlapping rings" begin
        for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
            @testset "$Criterion" test_static_nearly_overlapping_rings(Criterion)
        end
    end

    @testset "Static periodic self-reconnection" begin
        @testset "Nearly overlapping ellipses" begin
            for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
                @testset "$Criterion" test_static_periodic_nearly_overlapping_ellipses(Criterion)
            end
        end

        @testset "Overlapping ellipses" begin
            for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
                @testset "$Criterion" test_static_periodic_overlapping_ellipses(Criterion)
            end
        end

        @testset "Overlapping circles" begin
            for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
                @testset "$Criterion" test_static_periodic_overlapping_circles(Criterion)
            end
        end
    end

    @testset "Static: infinite line + ring" begin
        for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
            @testset "$Criterion" test_static_periodic_line_with_ring(Criterion)
        end
    end

    @testset "Static: two infinite lines" begin
        for Criterion in (ReconnectBasedOnDistance, ReconnectFast)
            @testset "$Criterion" test_static_periodic_line_with_line(Criterion)
        end
    end
end

# This can be useful for playing around with the tests.
if @isdefined(Makie)
    fig = Figure()
    Ls = 2π .* (1, 1, 1)
    # ax = Axis3(fig[1, 1]; aspect = :data)
    ax = Axis(fig[1, 1]; aspect = DataAspect())
    wireframe!(ax, Rect(0, 0, 0, Ls...); color = :grey, linewidth = 0.5)
    for f ∈ fs_orig
        p = plot!(ax, f; refinement = 8, linestyle = :dash, color = (:grey, 0.8), markersize = 8)
        inds = (firstindex(f) - 0):(lastindex(f) + 0)
        for i ∈ inds
            align = i ∈ eachindex(f) ? (:left, :bottom) : (:right, :top)
            text!(ax, f[i]; text = string(i), fontsize = 12, color = p.color, align)
        end
    end
    for f ∈ fs
        p = plot!(ax, f; refinement = 1, linestyle = :solid, markersize = 0, linewidth = 2)
        scatterlines!(ax, nodes(f).data; color = p.color, markersize = 8)
        inds = (firstindex(f) - 0):(lastindex(f) + 0)
        # for i ∈ inds
        #     align = i ∈ eachindex(f) ? (:left, :bottom) : (:right, :top)
        #     text!(ax, f[i]; text = string(i), fontsize = 12, color = p.color, align)
        # end
    end
    fig
end
