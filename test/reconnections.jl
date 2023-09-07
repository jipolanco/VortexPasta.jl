using Test
using VortexPasta.PredefinedCurves:
    define_curve, Ring, TrefoilKnot, Lemniscate, PeriodicLine
using VortexPasta.Filaments
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
using VortexPasta.Timestepping: VortexFilamentSolver
using StaticArrays
using Rotations
using UnicodePlots: lineplot
using LinearAlgebra: norm, I, ⋅
using JET: @test_opt
using FINUFFT: FINUFFT  # for JET only

# Note: this is total energy within a unit cell.
# One should normalise by the cell volume to get energy per unit mass [L²T⁻²].
function kinetic_energy_from_streamfunction(
        ψs_all::AbstractVector, fs::AbstractVector{<:AbstractFilament},
        Γ::Real;
        quad = nothing,
    )
    prefactor = Γ / 2
    E = zero(prefactor)
    interpolate_tangent_component_only = true
    for (f, ψs) ∈ zip(fs, ψs_all)
        ts = knots(f)
        if quad === nothing
            for i ∈ eachindex(segments(f))
                ψ⃗ = ψs[i]
                s⃗′ = f[i, Derivative(1)]
                δt = (ts[i + 1] - ts[i - 1]) / 2
                E += (ψ⃗ ⋅ s⃗′) * δt
            end
        else
            Xoff = Filaments.end_to_end_offset(f)
            ψ_int = Filaments.change_offset(similar(f), zero(Xoff))
            if interpolate_tangent_component_only
                for i ∈ eachindex(ψs)
                    ψ_t = ψs[i] ⋅ f[i, UnitTangent()]
                    ψ_int[i] = (ψ_t, 0, 0)  # we only care about the first component
                end
            else
                copy!(nodes(ψ_int), ψs)
            end
            update_coefficients!(ψ_int; knots = knots(f))
            E += integrate(f, quad) do f, i, ζ
                ψ⃗ = ψ_int(i, ζ)
                if interpolate_tangent_component_only
                    ψ⃗[1]
                else
                    s⃗′ = f(i, ζ, Derivative(1))
                    ψ⃗ ⋅ s⃗′
                end
            end
        end
    end
    E * prefactor
end

function kinetic_energy_from_streamfunction(iter::VortexFilamentSolver)
    (; fs, ψs,) = iter
    (; Γ,) = iter.prob.p.common
    quad = GaussLegendre(4)  # this doesn't seem to change much the results...
    kinetic_energy_from_streamfunction(ψs, fs, Γ; quad)
end

# Create a curve which resembles an 8 (or ∞).
# This specific curve is called the lemniscate of Bernouilli (https://en.wikipedia.org/wiki/Lemniscate_of_Bernoulli)
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

trefoil_scheme_dt(::Any) = trefoil_scheme_dt(RK4())  # default
trefoil_scheme_dt(::RK4) = 1

function test_trefoil_knot_reconnection(scheme = RK4())
    S = define_curve(TrefoilKnot(); translate = π, scale = π / 4)
    N = 96
    f = Filaments.init(S, ClosedFilament, N, CubicSplineMethod())
    fs_init = [f]
    l_min = minimum_knot_increment(fs_init)

    params_bs = let
        Ls = (1, 1, 1) .* 2π
        Ns = (1, 1, 1) .* 32
        kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
        α = kmax / 5
        rcut = 4 / α
        ParamsBiotSavart(;
            Γ = 2.0, α, a = 1e-6, Δ = 1/4, rcut, Ls, Ns,
            backend_short = NaiveShortRangeBackend(),
            backend_long = FINUFFTBackend(),
            quadrature_short = GaussLegendre(4),
            quadrature_long = GaussLegendre(4),
        )
    end

    tspan = (0.0, 1.0)  # ignored
    @test_opt VortexFilamentProblem(fs_init, tspan, params_bs)
    prob = @inferred VortexFilamentProblem(fs_init, tspan, params_bs)
    reconnect = ReconnectBasedOnDistance(1.5 * l_min)
    iter = @inferred init(
        prob, scheme;
        dt = 0.005 * trefoil_scheme_dt(scheme),
        refinement = RefineBasedOnSegmentLength(0.75 * l_min),
        # refinement = RefineBasedOnCurvature(0.4; ℓ_max = 1.5 * l_min, ℓ_min = 0.4 * l_min),
        reconnect,
        adaptivity = NoAdaptivity(),
    )
    (; fs,) = iter

    times = [iter.time.t]
    energy = [@inferred kinetic_energy_from_streamfunction(iter)]

    # TODO check that energy decreases? (currently, it slightly increases before the
    # reconnection...)
    Nf_prev = length(fs)
    n_reconnect = 0
    t_reconnect = 0.0
    @time for n = 2:500  # corresponds to the index in `times` array
        Timestepping.step!(iter)
        (; t, dt,) = iter.time
        push!(times, t)
        push!(energy, kinetic_energy_from_streamfunction(iter))
        Nf = length(fs)

        # let s = nameof(typeof(scheme))
        #     write_vtkhdf("trefoil_$(s)_$(n).hdf", fs; refinement = 8) do io
        #         io["Velocity"] = iter.vs
        #     end
        # end

        # First self-reconnection.
        if n_reconnect == 0 && Nf == 2
            t_reconnect = t
            n_reconnect = n
        end

        Nf_prev = Nf

        t > 1.95 && break
    end

    @test_opt ignored_modules=(Base, FINUFFT) step!(iter)

    @show t_reconnect
    @test n_reconnect > 0
    @test 1.65 < t_reconnect < 1.75  # for now this depends on a lot parameters...
    @test last(energy) < first(energy)

    let
        plt = lineplot(
            times, energy;
            xlabel = "Time", ylabel = "Energy", title = "Trefoil knot / $scheme",
        )
        println(plt)
    end

    nothing
end

@testset "Reconnections" begin
    @testset "Static self-reconnection: figure-eight knot" begin
        Az = 0.001
        curve = figure_eight_curve(; a = π / 4, Az, translate = (0, 0, 0))
        (; S, tlims,) = curve
        N = 32
        ζs = range(tlims...; length = 2N + 1)[2:2:2N]
        f = Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())

        l_min = minimum_knot_increment(f)

        @testset "Filaments.split!" begin
            i = length(f) ÷ 4
            j = 3i
            @test_opt ignored_modules=(Base,) Filaments.split!(copy(f), i, j)
            f1, f2 = @inferred Filaments.split!(copy(f), i, j)
            update_coefficients!.((f1, f2))
            @test length(f1) == length(f2) == length(f) ÷ 2
            @test f1 == f[i + 1:j]
            @test f2 == vcat(f[j + 1:end], f[begin:i])
        end

        @testset "reconnect_self!" begin
            crit = @inferred ReconnectBasedOnDistance(l_min / 2)
            fc = copy(f)
            fs_all = [fc]
            @test_opt ignored_modules=(Base,) Filaments.reconnect_self!(crit, fc, fs_all)
            f1 = Filaments.reconnect_self!(crit, fc, fs_all)
            @test f1 !== nothing       # reconnection happened
            @test length(fs_all) == 2  # filament split into two!
            fs_all[1] = f1  # one usually wants to replace the old `fc` filament, no longer valid
            f2 = fs_all[2]

            # Check that the reconnection happened at the right place
            i = length(f) ÷ 4
            j = 3i
            @test length(f1) == length(f2) == length(f) ÷ 2
            @test f1 == f[i + 1:j]
            @test f2 == vcat(f[j + 1:end], f[begin:i])
        end
    end

    @testset "Dynamic: trefoil knot" begin
        schemes = (RK4(), SanduMRI33a(8))
        @testset "Scheme: $scheme" for scheme ∈ schemes
            test_trefoil_knot_reconnection(scheme)
        end
    end

    @testset "Static: antiparallel rings" begin
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
        crit = ReconnectBasedOnDistance(1.2 * d_min_nodes)

        @testset "Filaments.reconnect_other!" begin
            f, g = copy.(fs_orig)
            h = @inferred Nothing Filaments.reconnect_other!(crit, f, g)
            @test h !== nothing  # reconnection happened!
            @test length(h) == sum(length, fs_orig)
            @test maximum_knot_increment(h) < 2 * l_min  # check that there are no crazy jumps
        end

        @testset "Filaments.reconnect!" begin
            fs = collect(copy.(fs_orig))
            nrec = @inferred Filaments.reconnect!(crit, fs)
            @test nrec == 1  # one reconnection
            @test length(fs) == 1
            @test length(fs[1]) == sum(length, fs_orig)
            @test maximum_knot_increment(fs) < 2 * l_min  # check that there are no crazy jumps
        end
    end

    @testset "Static: nearly overlapping rings" begin
        # Test two nearly overlapping rings (with some small offset in z) which should reconnect at two points.
        # We should end up with two separate vortices.
        # Internally, this happens in two steps:
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
        crit = ReconnectBasedOnDistance(1.2 * d_min_nodes)

        fs = collect(copy.(fs_orig))
        nrec = @inferred Filaments.reconnect!(crit, fs)
        @test nrec == 2
        @test length(fs) == 2
        @test fs[1] != fs_orig[1]  # reconnection happened
        @test fs[2] != fs_orig[2]  # reconnection happened
        @test length(fs[1]) ≠ length(fs[2])  # the resulting vortices have different sizes
    end

    @testset "Static periodic self-reconnection" begin
        @testset "Nearly overlapping ellipses" begin
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

            fs = [copy(f_orig)]
            crit = ReconnectBasedOnDistance(4 * ϵ)
            f1 = Filaments.reconnect_self!(crit, first(fs), fs; periods)
            @test f1 !== nothing
            @test length(fs) == 2  # filament reconnected onto 2
            fs[1] = f1  # replace now invalid filament with new filament
            @test abs(end_to_end_offset(fs[1])[1]) ≈ 2π  # infinite line in x
            @test abs(end_to_end_offset(fs[2])[1]) ≈ 2π  # infinite line in x
            foreach(Filaments.update_coefficients!, fs)
            @test all(f -> no_jumps(f, 2 * l_min), fs)
        end

        @testset "Overlapping ellipses" begin
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
            crit = ReconnectBasedOnDistance(l_min / 2)
            nrec = @inferred Filaments.reconnect!(crit, fs; periods)
            @test nrec == 2
            @test length(fs) == 3
            @test sum(length, fs) == N  # number of nodes didn't change

            # Note that the order of the filaments is arbitrary and implementation-dependent...
            @test sum(f -> norm(end_to_end_offset(f)), fs) ≈ 4π  # two infinite lines
            @test end_to_end_offset(fs[1]) ≈ Vec3(-2π, 0, 0)  # infinite line
            @test end_to_end_offset(fs[2]) == Vec3(0, 0, 0)   # small loop
            @test end_to_end_offset(fs[3]) ≈ Vec3(+2π, 0, 0)  # infinite line
            @test all(f -> no_jumps(f, 2.1 * l_min), fs)
        end

        @testset "Overlapping circles" begin
            # Slightly more complex case of overlapping periodic circles.
            # Each circle is reconnected into 3 closed curves.
            # Right now, this requires two passes of Filaments.reconnect! (shouldn't really
            # matter in practice...).
            periods = 2π .* (1, 1, 1)
            R = π * 1.2
            N = 64
            ring = ring_curve(; scale = R, translate = π * Vec3(1, 1, 1),)
            f_orig = let
                (; S, tlims,) = ring
                ζs = range(tlims...; length = 2N + 1)[2:2:2N]
                Filaments.init(ClosedFilament, S.(ζs), FiniteDiffMethod(2))
            end

            l_min = minimum_knot_increment(f_orig)
            fs = [copy(f_orig)]
            crit = ReconnectBasedOnDistance(l_min / 2)

            # First pass: the ring splits into 2 infinite lines and a closed curve (very
            # similar to the overlapping ellipse case).
            nrec = Filaments.reconnect!(crit, fs; periods)
            @test nrec == 2
            @test length(fs) == 3
            @test sum(length, fs) == N  # number of nodes didn't change
            @test sum(f -> norm(end_to_end_offset(f)), fs) ≈ 4π  # two infinite lines
            @test all(f -> no_jumps(f, 2.1 * l_min), fs)

            # Second pass: the two infinite lines merge and then split, forming two closed lines.
            nrec = Filaments.reconnect!(crit, fs; periods)
            @test nrec == 2
            @test length(fs) == 3
            @test sum(length, fs) == N  # number of nodes didn't change
            @test sum(f -> norm(end_to_end_offset(f)), fs) < 1e-12  # no infinite lines
            @test all(f -> no_jumps(f, 2.1 * l_min), fs)
        end
    end

    @testset "Static: infinite line + ring" begin
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
        crit = ReconnectBasedOnDistance(0.5 * l_min)
        nrec = Filaments.reconnect!(crit, fs)
        @test nrec == 2
        for f ∈ fs
            Filaments.fold_periodic!(f, periods)
            Filaments.update_coefficients!(f)
        end
        @test length(fs) == 2
        @test fs == fs_manual  # compare with manually merged vortices
    end

    @testset "Static: two infinite lines" begin
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
        crit = ReconnectBasedOnDistance(l_min)
        periods = 2π .* (1, 1, 1)
        fs = copy.(fs_orig)
        nrec = Filaments.reconnect!(crit, fs; periods)
        @test nrec == 1
        @test length(fs) == 1
        @test fs == fs_manual  # compare with manually merged vortices
    end

end

# This can be useful for playing around with the tests.
if false && @isdefined(Makie)
    fig = Figure()
    Ls = 2π .* (1, 1, 1)
    ax = Axis3(fig[1, 1]; aspect = :data)
    wireframe!(ax, Rect(0, 0, 0, Ls...); color = :grey, linewidth = 0.5)
    # for f ∈ fs_orig
    #     plot!(ax, f; refinement = 8, linestyle = :dash, color = (:grey, 0.5), markersize = 8)
    # end
    for f ∈ fs
        p = plot!(ax, f; refinement = 8, linestyle = :solid)
        scatter!(ax, nodes(f).data; color = p.color)
        inds = (firstindex(f) - 1):(lastindex(f) + 1)
        for i ∈ inds
            align = i ∈ eachindex(f) ? (:left, :bottom) : (:right, :top)
            text!(ax, f[i]; text = string(i), fontsize = 16, color = p.color, align)
        end
    end
    fig
end
