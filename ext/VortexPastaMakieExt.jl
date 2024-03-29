module VortexPastaMakieExt

using VortexPasta.Filaments

using MakieCore: MakieCore
using Observables: Observable, @map
using LinearAlgebra: normalize

const MaybeObservable{T} = Union{Observable{T}, T}

MakieCore.plot(f::MaybeObservable{<:AbstractFilament}, args...; kws...) =
    filamentplot(f, args...; kws...)
MakieCore.plot!(ax, f::MaybeObservable{<:AbstractFilament}, args...; kws...) =
    filamentplot!(ax, f, args...; kws...)

# This macro defines the functions `filamentplot` and `filamentplot!` and the type FilamentPlot,
# among other things.
MakieCore.@recipe(FilamentPlot) do scene
    MakieCore.Attributes(
        refinement = 1,
        marker = :circle,
        markersize = 10.0f0,
        linewidth = 1.5f0,
        linestyle = :solid,
        cycle = [:color],  # this gives the default colour cycle (see e.g. docs for Makie.lines)
        color = :black,
        colormap = :viridis,
        colorrange = MakieCore.Automatic(),
        markercolor = nothing,  # defaults to `color`
        tangents = false,
        curvatures = false,
        tangentcolor = nothing,
        curvaturecolor = nothing,
        vectorpos = 0.0,
        velocitycolor = nothing,
        periods = (nothing, nothing, nothing),

        # Arrow attributes
        arrowscale = 1.0f0,
        arrowwidth = MakieCore.Automatic(),
        arrowsize = MakieCore.Automatic(),
    )
end

const PlotFilamentOnly = FilamentPlot{<:Tuple{<:AbstractFilament}}
const PlotFilamentAndVelocities = FilamentPlot{<:Tuple{<:AbstractFilament, <:AbstractVector{<:Vec3}}}

# For some reason the second argument (`len`) is also needed in dispatch (not really mentioned in the Makie docs...).
MakieCore.argument_names(::Type{<:PlotFilamentOnly}, len::Integer) = (:filament,)
MakieCore.argument_names(::Type{<:PlotFilamentAndVelocities}, len::Integer) = (:filament, :velocities)

fixfirst(f::F, x) where {F} = (ys...) -> f(x, ys...)

# We use the fact that `Makie.lift` (not included in MakieCore) is the same as `map`.
function MakieCore.plot!(p::FilamentPlot)
    argnames = MakieCore.argument_names(p)
    f = p.filament :: Observable{<:AbstractFilament}
    v = (:velocities ∈ argnames) ? p.velocities : nothing
    Xs_line = @map _refine_filament_for_plotting(&f, &p.refinement, &p.periods)
    MakieCore.lines!(
        p, Xs_line;
        color = p.color, linewidth = p.linewidth, linestyle = p.linestyle,
        colormap = p.colormap, colorrange = p.colorrange,
    )
    Xs_nodes = map(f, p.periods) do f, Ls
        if all(isnothing, Ls)
            nodes(f)
        else
            map(s⃗ -> to_main_periodic_cell(s⃗, Ls), nodes(f))
        end
    end
    MakieCore.scatter!(
        p, Xs_nodes;
        color = @map(something(&p.markercolor, &p.color)),
        marker = p.marker, markersize = p.markersize,
        colormap = p.colormap, colorrange = p.colorrange,
    )
    arrow_attrs = @map (&p.arrowscale, &p.arrowwidth, &p.arrowsize)  # new observable from the 3 observables
    let
        tangentcolor = @map something(&p.tangentcolor, &p.color)
        @map _plot_tangents!(
            p, f, &p.tangents, &p.vectorpos, &p.periods, &arrow_attrs;
            color = &tangentcolor, colormap = &p.colormap, colorrange = &p.colorrange,
        )
    end
    let
        curvaturecolor = @map something(&p.curvaturecolor, &p.color)
        @map _plot_curvatures!(
            p, f, &p.curvatures, &p.vectorpos, &p.periods, &arrow_attrs;
            color = &curvaturecolor, colormap = &p.colormap, colorrange = &p.colorrange,
        )
    end
    if v !== nothing
        velocitycolor = @map something(&p.velocitycolor, &p.color)
        @map _plot_velocities!(
            p, f, v, &p.periods, &arrow_attrs;
            color = &velocitycolor, colormap = &p.colormap, colorrange = &p.colorrange,
        )
    end
    p
end

to_main_periodic_cell(x⃗, Ls::Tuple) = oftype(x⃗, map(to_main_periodic_cell, x⃗, Ls))
to_main_periodic_cell(x, L::Nothing) = x  # don't do anything

function to_main_periodic_cell(x, L::Real)
    while x > L
        x -= L
    end
    while x < 0
        x += L
    end
    x
end

is_jump(x⃗, x⃗_prev, Lhs::Tuple) = any(splat(is_jump), zip(x⃗, x⃗_prev, Lhs))
is_jump(x, x_prev, Lh::Nothing) = false
is_jump(x, x_prev, Lh::Real) = abs(x - x_prev) ≥ Lh

function _refine_filament_for_plotting(
        f::ClosedFilament, refinement::Int, Ls::Tuple,
    )
    Xs_nodes = nodes(f)
    refinement ≥ 1 || error("refinement must be ≥ 1")
    # Note: if `Ls` contains values different from `nothing`, then lines may need to
    # be "broken" by inserting NaNs, and we don't know a priori how many points we need to
    # return.
    with_periods = any(!isnothing, Ls)
    T = eltype(Xs_nodes)
    x⃗_nan = fill(NaN, T)::T
    Xs = similar(Xs_nodes, 0)
    N = refinement * length(f) + 1  # the +1 is to close the loop
    sizehint!(Xs, N)  # this length estimate is only true when we don't need to insert NaNs
    n = 0
    subinds = range(0, 1; length = refinement + 1)[1:refinement]
    Lhs = map(L -> L === nothing ? nothing : L / 2, Ls)  # half periods
    x⃗_prev = to_main_periodic_cell(first(Xs_nodes), Ls)
    njumps = 0
    for i ∈ eachindex(f)
        for ζ ∈ subinds
            x⃗ = to_main_periodic_cell(f(i, ζ), Ls)
            # Check if we "jumped" from the previous position.
            # Insert a NaN point in that case.
            if with_periods
                if is_jump(x⃗, x⃗_prev, Lhs)
                    njumps += 1
                    push!(Xs, x⃗_nan)
                end
            end
            push!(Xs, x⃗)
            x⃗_prev = x⃗
        end
    end
    @assert length(Xs) == refinement * length(f) + njumps
    # Close the loop.
    let x⃗ = to_main_periodic_cell(f(lastindex(f), 1.0), Ls)
        if !is_jump(x⃗, x⃗_prev, Lhs)
            push!(Xs, x⃗)
        end
    end
    Xs
end

function _tangents_for_arrows(f::AbstractFilament, vectorpos::Real, Ls)
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = f(i, vectorpos, UnitTangent())
        for n ∈ eachindex(x)
            Xs[n][i] = to_main_periodic_cell(x[n], Ls[n])
            Vs[n][i] = v[n]
        end
    end
    (Xs..., Vs...)
end

function _plot_arrows!(p, arrow_attrs, args...; kwargs...)
    MakieCore.arrows!(p, args...; kwargs..., _arrow_kwargs(arrow_attrs)...)
end

function _plot_tangents!(
        p::FilamentPlot, f::Observable{<:AbstractFilament}, tangents::Bool,
        vectorpos::Real, periods::Tuple, arrow_attrs::Tuple;
        kwargs...,
    )
    tangents || return p
    data = @map _tangents_for_arrows(&f, vectorpos, periods)
    args = ntuple(i -> @map((&data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrow_attrs, args...; kwargs...)
    p
end

function _curvatures_for_arrows(f::AbstractFilament, vectorpos::Real, Ls)
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = f(i, vectorpos, CurvatureVector())
        # ρ = norm(v)  # curvature (= 1 / R)
        for n ∈ eachindex(x)
            Xs[n][i] = to_main_periodic_cell(x[n], Ls[n])
            Vs[n][i] = v[n]
        end
    end
    (Xs..., Vs...)
end

function _plot_curvatures!(
        p::FilamentPlot, f::Observable{<:AbstractFilament}, curvatures::Bool,
        vectorpos::Real, periods::Tuple, arrow_attrs::Tuple;
        kwargs...,
    )
    curvatures || return p
    data = @map _curvatures_for_arrows(&f, vectorpos, periods)
    args = ntuple(i -> @map((&data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrow_attrs, args...; kwargs...)
    p
end

function _velocities_for_arrows(f::AbstractFilament, v::AbstractVector, Ls)
    length(v) == length(f) || throw(DimensionMismatch("wrong length of velocity vector"))
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f), n ∈ eachindex(Xs)
        Xs[n][i] = to_main_periodic_cell(f[i][n], Ls[n])
        Vs[n][i] = v[i][n]
    end
    (Xs..., Vs...)
end

function _plot_velocities!(
        p::FilamentPlot, f::Observable{<:AbstractFilament},
        v::Observable{<:AbstractVector}, periods::Tuple, arrow_attrs::Tuple;
        kwargs...,
    )
    data = @map _velocities_for_arrows(&f, &v, periods)
    args = ntuple(i -> @map((&data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrow_attrs, args...; kwargs...)
    p
end

# Convert input arrow arguments to arguments for MakieCore.arrows!.
function _arrow_kwargs(arrow_attrs::Tuple)
    arrowscale, arrowwidth, arrowsize = arrow_attrs
    (;
        arrowsize,
        linewidth = arrowwidth,
        lengthscale = arrowscale,
    )
end

end
