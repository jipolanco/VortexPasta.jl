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

- `color = :black`

- `linewidth = 1.5f0`

- `linestyle = :solid`

- `markercolor = nothing` (`nothing` → same as `color`)

- `marker = :circle`

- `markersize = 10.0f0`

- `colormap = :viridis`

- `colorrange = MakieCore.Automatic()`

### Arrow arguments

The following are used when plotting arrows (tangents, curvatures, velocities, …):

- `arrowscale = 1.0f0` allows to scale vectors (controls their *length*).
  Corresponds to `lengthscale` in `Makie.arrows`;

- `arrowsize = MakieCore.Automatic()` controls the head size.
  It has the same name in `Makie.arrows`.
  It should be a tuple `(sx, sy, sz)`, where the first 2 set the width of the cone, and `sz`
  its height;

- `arrowwidth = MakieCore.Automatic()` sets the linewidth of the arrow.
  Corresponds to `linewidth` in `Makie.arrows`.

See also the [Makie docs on arrows](https://docs.makie.org/stable/reference/plots/arrows/index.html).

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
