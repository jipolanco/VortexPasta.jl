module VortexPastaMakieExt

using VortexPasta.Filaments

using MakieCore: MakieCore
using Observables: Observable, @map
using LinearAlgebra: normalize

MakieCore.plot(f::AbstractFilament, args...; kws...) = filamentplot(f, args...; kws...)
MakieCore.plot!(ax, f::AbstractFilament, args...; kws...) = filamentplot!(ax, f, args...; kws...)

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
        markercolor = nothing,  # defaults to `color`
        tangents = false,
        curvatures = false,
        tangentcolor = nothing,
        curvaturecolor = nothing,
        vectorpos = 0.0,
        velocitycolor = nothing,

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
    Xs_nodes = @map _select_points_to_plot(&f)
    Xs_line = @map Filaments._refine_filament(&f, &p.refinement)
    MakieCore.lines!(
        p, Xs_line;
        color = p.color, linewidth = p.linewidth, linestyle = p.linestyle,
        colormap = p.colormap,
    )
    MakieCore.scatter!(
        p, Xs_nodes;
        color = @map(something(&p.markercolor, &p.color)),
        marker = p.marker, markersize = p.markersize,
        colormap = p.colormap,
    )
    arrow_attrs = @map (&p.arrowscale, &p.arrowwidth, &p.arrowsize)  # new observable from the 3 observables
    let
        tangentcolor = @map something(&p.tangentcolor, &p.color)
        @map _plot_tangents!(
            p, f, &p.tangents, &p.vectorpos, &tangentcolor, &p.colormap, &arrow_attrs,
        )
    end
    let
        curvaturecolor = @map something(&p.curvaturecolor, &p.color)
        @map _plot_curvatures!(
            p, f, &p.curvatures, &p.vectorpos, &curvaturecolor, &p.colormap, &arrow_attrs,
        )
    end
    if v !== nothing
        velocitycolor = @map something(&p.velocitycolor, &p.color)
        @map _plot_velocities!(p, f, v, &velocitycolor, &p.colormap, &arrow_attrs)
    end
    p
end

# Make sure we include `end + 1` point to close the loop.
_select_points_to_plot(f::ClosedFilament) = nodes(f)[begin:end + 1]

function _tangents_for_arrows(f::AbstractFilament, vectorpos::Real)
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = f(i, vectorpos, UnitTangent())
        for n ∈ eachindex(x)
            Xs[n][i] = x[n]
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
        vectorpos::Real, color, colormap, arrow_attrs::Tuple,
    )
    tangents || return p
    data = @map _tangents_for_arrows(&f, vectorpos)
    args = ntuple(i -> @map((&data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrow_attrs, args...; color, colormap)
    p
end

function _curvatures_for_arrows(f::AbstractFilament, vectorpos::Real)
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = f(i, vectorpos, CurvatureVector())
        # ρ = norm(v)  # curvature (= 1 / R)
        for n ∈ eachindex(x)
            Xs[n][i] = x[n]
            Vs[n][i] = v[n]
        end
    end
    (Xs..., Vs...)
end

function _plot_curvatures!(
        p::FilamentPlot, f::Observable{<:AbstractFilament}, curvatures::Bool,
        vectorpos::Real,
        color,
        colormap,
        arrow_attrs::Tuple,
    )
    curvatures || return p
    data = @map _curvatures_for_arrows(&f, vectorpos)
    args = ntuple(i -> @map((&data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrow_attrs, args...; color, colormap)
    p
end

function _velocities_for_arrows(f::AbstractFilament, v::AbstractVector)
    length(v) == length(f) || throw(DimensionMismatch("wrong length of velocity vector"))
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f), n ∈ eachindex(Xs)
        Xs[n][i] = f[i][n]
        Vs[n][i] = v[i][n]
    end
    (Xs..., Vs...)
end

function _plot_velocities!(
        p::FilamentPlot, f::Observable{<:AbstractFilament},
        v::Observable{<:AbstractVector}, color, colormap,
        arrow_attrs::Tuple,
    )
    data = @map _velocities_for_arrows(&f, &v)
    args = ntuple(i -> @map((&data)[i]), Val(6))  # convert Observable of tuples to tuple of Observables
    _plot_arrows!(p, arrow_attrs, args...; color, colormap)
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
