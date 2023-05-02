using MakieCore: MakieCore
using Observables: Observable, @map
using LinearAlgebra: normalize

export filamentplot, filamentplot!

"""
    filamentplot(f::AbstractFilament, [velocities]; kws...)
    MakieCore.plot(f::AbstractFilament, [velocities]; kws...)

Plot a filament using Makie.jl.

Example usage:

```julia
using GLMakie
plot(f; refinement = 4)  # f is a filament
```

See [`filamentplot!`](@ref) for details and for optional keyword arguments.
"""
function filamentplot end

"""
    filamentplot!([ax,] f::AbstractFilament, [velocities]; kws...)
    MakieCore.plot!(ax, f::AbstractFilament, [velocities]; kws...)

Plot filament onto existent 3D axis.

The first argument should typically be an `Axis3` or an `LScene`.

Example usage:

```julia
using GLMakie
fig = Figure()
ax = Axis3(fig[1, 1])
plot!(ax, f)  # f is a filament
```

## Optional arguments and their defaults

- `refinement::Int = 1`: level of refinement of the curves (must be ≥ 1)

- `color = :black`

- `linewidth = 1.5f0`

- `markercolor = nothing` (`nothing` → same as `color`)

- `marker = :circle`

- `markersize = 10.0f0`

- `arrowscale = 1.0f0` allows to scale vectors (tangents, curvatures, velocities, …)

- `colormap = :viridis`

### Plotting tangent and curvature vectors

Tangent and curvature vectors can be optionally plotted via the `tangents` and
`curvatures` arguments. A single vector will be plotted for each filament segments.
By default, vectors are evaluated at filament nodes, but one can also evaluate
them in-between nodes using the `vectorpos` argument.

- `tangents::Bool = false`: plot unit tangent vectors.

- `curvatures::Bool = false`: plot curvature vectors. Note that the magnitude
  is the local curvature ``ρ = 1 / R``, where ``R`` is the curvature *radius*.

- `tangentcolor = nothing`

- `curvaturecolor = nothing`

- `vectorpos = 0.0`: relative vector positions within each segment. Must be in ``[0, 1]``.

### Plotting velocities of filament nodes

Similarly, it is possible to plot vector quantities attached to filament nodes,
such as filament velocities. For this pass a vector of velocities as a
positional argument after the filament `f`.

Associated keyword arguments:

- `velocitycolor = nothing` colour of velocity vectors

"""
function filamentplot! end

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
        cycle = [:color],  # this gives the default colour cycle (see e.g. docs for Makie.lines)
        color = :black,
        colormap = :viridis,
        markercolor = nothing,  # defaults to `color`
        arrowscale = 1.0f0,
        tangents = false,
        curvatures = false,
        tangentcolor = nothing,
        curvaturecolor = nothing,
        vectorpos = 0.0,
        velocitycolor = nothing,
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
    Xs_line = @map _refine_filament(&f, &p.refinement)
    MakieCore.lines!(
        p, Xs_line;
        color = p.color, linewidth = p.linewidth,
        colormap = p.colormap,
    )
    MakieCore.scatter!(
        p, Xs_nodes;
        color = @map(something(&p.markercolor, &p.color)),
        marker = p.marker, markersize = p.markersize,
        colormap = p.colormap,
    )
    arrowscale = p.arrowscale
    let
        tangentcolor = @map something(&p.tangentcolor, &p.color)
        @map _plot_tangents!(
            p, &f, &p.tangents, &p.vectorpos, &tangentcolor, &p.colormap, &arrowscale,
        )
    end
    let
        curvaturecolor = @map something(&p.curvaturecolor, &p.color)
        @map _plot_curvatures!(
            p, &f, &p.curvatures, &p.vectorpos, &curvaturecolor, &p.colormap, &arrowscale,
        )
    end
    if v !== nothing
        velocitycolor = @map something(&p.velocitycolor, &p.color)
        @map _plot_velocities!(p, &f, &v, &velocitycolor, &p.colormap, &arrowscale)
    end
    p
end

# Make sure we include `end + 1` point to close the loop.
_select_points_to_plot(f::ClosedFilament) = nodes(f)[begin:end + 1]

function _plot_tangents!(
        p::FilamentPlot, f::ClosedFilament, tangents::Bool,
        vectorpos::Real, tangentcolor, colormap, arrowscale::Real,
    )
    tangents || return p
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        v = normalize(f(i, vectorpos, Derivative(1)))
        for n ∈ eachindex(x)
            Xs[n][i] = x[n]
            Vs[n][i] = v[n]
        end
    end
    MakieCore.arrows!(
        p, Xs..., Vs...;
        color = tangentcolor,
        colormap,
        _arrow_kwargs(; arrowscale,)...,
    )
    p
end

function _plot_curvatures!(
        p::FilamentPlot, f::ClosedFilament, curvatures::Bool,
        vectorpos::Real,
        curvaturecolor,
        colormap,
        arrowscale::Real,
    )
    curvatures || return p
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f)
        x = f(i, vectorpos, Derivative(0))
        x′ = f(i, vectorpos, Derivative(1))
        x″ = f(i, vectorpos, Derivative(2))
        _, v = normalise_derivatives(x′, x″)
        # ρ = norm(v)  # curvature (= 1 / R)
        for n ∈ eachindex(x)
            Xs[n][i] = x[n]
            Vs[n][i] = v[n]
        end
    end
    MakieCore.arrows!(
        p, Xs..., Vs...;
        color = curvaturecolor,
        colormap,
        _arrow_kwargs(; arrowscale,)...,
    )
    p
end

function _plot_velocities!(
        p::FilamentPlot, f::AbstractFilament, v::AbstractVector, velocitycolor,
        colormap,
        arrowscale::Real,
    )
    length(v) == length(f) || throw(DimensionMismatch("wront length of velocity vector"))
    N = length(f)
    Xs = ntuple(_ -> Vector{Float32}(undef, N), 3)  # layout accepted by Makie.arrows
    Vs = map(similar, Xs)
    for i ∈ eachindex(f), n ∈ eachindex(Xs)
        Xs[n][i] = f[i][n]
        Vs[n][i] = v[i][n]
    end
    MakieCore.arrows!(
        p, Xs..., Vs...;
        color = velocitycolor,
        colormap,
        _arrow_kwargs(; arrowscale,)...,
    )
    p
end

_arrow_kwargs(; arrowscale,) = (;
    arrowsize = (0.15, 0.15, 0.25),
    linewidth = 0.04,
    lengthscale = arrowscale,
)

function _refine_filament(f::ClosedFilament, refinement::Int)
    Xs_nodes = nodes(f)
    refinement ≥ 1 || error("refinement must be ≥ 1")
    refinement == 1 && return Xs_nodes[begin:end + 1]
    N = refinement * length(f) + 1  # the +1 is to close the loop
    Xs = similar(Xs_nodes, N)
    n = 0
    subinds = range(0, 1; length = refinement + 1)[1:refinement]
    for i ∈ eachindex(f)
        for ζ ∈ subinds
            n += 1
            Xs[n] = f(i, ζ)
        end
    end
    @assert n == refinement * length(f)
    Xs[n + 1] = f(lastindex(f), 1.0)  # close the loop
    Xs
end
