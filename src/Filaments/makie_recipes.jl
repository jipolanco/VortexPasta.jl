using MakieCore: MakieCore
using LinearAlgebra: normalize

export filamentplot, filamentplot!

"""
    filamentplot(f::AbstractFilament; kws...)
    MakieCore.plot(f::AbstractFilament; kws...)

Plot a filament using Makie.jl.

Example usage:

```julia
using GLMakie
plot(f; axis = (type = Axis3,), refinement = 4)  # f is a filament
```

See [`filamentplot!`](@ref) for details and for optional keyword arguments.
"""
function filamentplot end

"""
    filamentplot!(ax, f::AbstractFilament; kws...)
    MakieCore.plot!(ax, f::AbstractFilament; kws...)

Plot filament onto 3D axis.

The first argument should ideally be an `Axis3`.

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
"""
function filamentplot! end

MakieCore.plot(f::AbstractFilament; kws...) = filamentplot(f; kws...)
MakieCore.plot!(ax, f::AbstractFilament; kws...) = filamentplot!(ax, f; kws...)

# This macro defines the functions `filamentplot` and `filamentplot!` and the type FilamentPlot,
# among other things.
MakieCore.@recipe(FilamentPlot, filament) do scene
    MakieCore.Attributes(
        refinement = 1,
        marker = :circle,
        markersize = 10.0f0,
        linewidth = 1.5f0,
        color = :black,
        markercolor = nothing,  # defaults to `color`
        tangents = false,
        curvatures = false,
        tangentcolor = nothing,
        curvaturecolor = nothing,
        vectorpos = 0.0,
    )
end

# We use the fact that `Makie.lift` (not included in MakieCore) is the same as `map`.
function MakieCore.plot!(p::FilamentPlot)
    f = p.filament
    f[] isa AbstractFilament || error("expected an AbstractFilament as first argument")
    Xs_nodes = map(_select_points_to_plot, f)
    Xs_line = map(_refine_filament, f, p.refinement)
    MakieCore.lines!(
        p, Xs_line;
        color = p.color, linewidth = p.linewidth,
    )
    MakieCore.scatter!(
        p, Xs_nodes;
        color = map(something, p.markercolor, p.color),
        marker = p.marker, markersize = p.markersize,
    )
    tangentcolor = map(something, p.tangentcolor, p.color)
    map(f, p.tangents, p.vectorpos, tangentcolor) do (args...)
        _plot_tangents!(p, args...)
    end
    curvaturecolor = map(something, p.curvaturecolor, p.color)
    map(f, p.curvatures, p.vectorpos, curvaturecolor) do (args...)
        _plot_curvatures!(p, args...)
    end
    p
end

# Make sure we include `end + 1` point to close the loop.
_select_points_to_plot(f::ClosedFilament) = points(f)[begin:end + 1]

function _plot_tangents!(
        p::FilamentPlot, f::ClosedFilament, tangents::Bool,
        vectorpos::Real,
        tangentcolor,
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
        _arrow_kwargs()...,
    )
    p
end

function _plot_curvatures!(
        p::FilamentPlot, f::ClosedFilament, curvatures::Bool,
        vectorpos::Real,
        curvaturecolor,
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
        _arrow_kwargs()...,
    )
    p
end

_arrow_kwargs() = (;
    arrowsize = (0.15, 0.15, 0.25),
    linewidth = 0.04,
)

function _refine_filament(f::ClosedFilament, refinement::Int)
    Xs_nodes = points(f)
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
