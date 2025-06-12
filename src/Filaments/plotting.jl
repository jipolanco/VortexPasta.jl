# These are actually defined as a package extension when Makie.jl is also installed.
# See ext/VortexPastaMakieExt.jl.
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

- `periods = (nothing, nothing, nothing)`. This can be a tuple of values representing the
  domain period, for instance `(2π, 2π, 2π)`. In that case, filaments will be "broken" when
  they exit the main unit cell, so that all vortex elements are within the cell.

- `color = :black`

- `linewidth = 1.5f0`

- `linestyle = :solid`

- `markercolor = nothing` (`nothing` → same as `color`)

- `marker = :circle`

- `markersize = 10.0f0`

- `colormap = :viridis`

- `colorrange = MakieCore.Automatic()`

### Arrow properties

One can pass an `arrows3d` named tuple to set arrow properties. For example:

`arrows3d = (; shaftlength = 0.6, lengthscale = 2.0,)`

See the [Makie docs](https://docs.makie.org/stable/reference/plots/arrows#Arrows3D) for a
full list of available options.

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
