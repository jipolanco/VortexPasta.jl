# # Getting started
#
# The following tutorial goes through the simulation of a **single vortex ring** propagating
# due to its self-induced propulsion.
#
# The aim is to show how one can:
#
# 1. Discretise a vortex filament using a spatial curve;
# 2. compute things like the self-induced vortex velocity;
# 3. simulate the motion of the vortex filament over a time period.
#
# We assume that VortexPasta and GLMakie have been installed following the instructions in
# the [Installation](@ref) section.
# The code in this tutorial should be executed in the same local environment where those
# packages were installed.
#
# ## Initialising a vortex ring
#
# The first thing to do is to define a circular vortex ring.
# In VortexPasta, curves are initialised using the [`Filaments.init`](@ref) function defined
# in the [`VortexPasta.Filaments`](@ref) submodule.
#
# The most straightforward way of defining an arbitrary curve is by first defining a set of
# discretisation points (which we also call *nodes*) and then passing that to `Filaments.init`.
# In VortexPasta, a point in 3D space is represented by the [`Vec3`](@ref
# VortexPasta.BasicTypes.Vec3) type (which is nothing else than an `SVector` from the
# [StaticArrays.jl](https://github.com/JuliaArrays/StaticArrays.jl) package).
#
# Let's define a set of points discretising a circular ring of radius ``R`` living on the
# plane ``z = 1`` and centred at ``\bm{x}_0 = [3, 3, 1]``:

using VortexPasta
using VortexPasta.Filaments

R = 0.5  # radius of the circular ring
N = 16   # number of discretisation points
x⃗₀ = Vec3(3.0, 3.0, 1.0)  # ring centre
θs = range(0, 2π; length = N + 1)[1:N]  # discretisation angles (we exclude θ = 2π)
points = [x⃗₀ + R * Vec3(cos(θ), sin(θ), 0) for θ ∈ θs]

# Now that we have defined a set of points, we can create a filament using `Filaments.init`.
# Note that this function requires choosing a **discretisation method** which will be used to
# interpolate the curve in-between nodes and to estimate curve derivatives.
# Here we use the [`CubicSplineMethod`](@ref), which represents curves as **periodic cubic
# splines**:

f = Filaments.init(ClosedFilament, points, CubicSplineMethod())

# Another possible discretisation option is the [`FiniteDiffMethod`](@ref), which estimates
# derivatives at discretisation points using **finite differences** (based on the locations of
# neighbouring points), and performs Hermite interpolations using those derivatives to
# reconstruct the curve in-between nodes.
#
# Besides, note that the first argument ([`ClosedFilament`](@ref)) is mandatory and is only
# there to make sure that we're creating a *closed* (as opposed to an open-ended) filament.
# This means in particular that the filament will be automatically "closed" when evaluated
# outside of its range of definition `1:N`.
# For example, one has `f[1] == f[N + 1]`:

@show f[1] f[N] f[N + 1]
nothing  # hide

# ### Getting geometric information out of a filament
#
# Above we have used the `f[i]` syntax to obtain the location of a discretisation point.
# As mentioned, curve locations can also be evaluated in-between discretisation points
# using **interpolation**.
# For instance, to evaluate the curve location at some point in-between nodes `i` and `i +
# 1`, one can call `f(i, ζ)`, where `ζ` is some real value between 0 and 1.
# In particular, `f(i, 0.5)` gives an estimation of the midpoint between both nodes,
# while `f(i, 0.0)` and `f(i, 1.0)` respectively correspond to `f[i]` and `f[i + 1]`:

i = 2
@show f[i] f(i, 0.0) f(i, 0.5) f(i, 1.0) f[i + 1]
nothing  # hide

# A similar syntax can be used to obtain **derivatives** at discretisation points or
# in-between them:

@show f[i] f[i, Derivative(1)] f[i, Derivative(2)]
@show f(i, 0.5) f(i, 0.5, Derivative(1)) f(i, 0.5, Derivative(2))
nothing  # hide

# These derivatives assume that curves are parametrised as ``\bm{s}(t)`` for ``t ∈ [0, T]``,
# and are computed with respect to this parameter.
# In practice, the parameter ``t`` roughly corresponds to the integrated arc length (and
# thus ``T`` is a rough estimation of the total curve length), but this should never be
# assumed.
# In particular, the first derivative ``∂\bm{s}/∂t`` is tangent to the curve but it's not
# necessarily unitary (which would be the case if ``t`` was the actual arc length ``ξ``).
#
# In some cases, one may want to directly obtain geometrically-relevant quantities such as
# the unit tangent or the curvature vector (see [Geometric quantities](@ref) for definitions
# and a list of possible quantities):

t̂ = f[i, UnitTangent()]       # `t̂` can be obtained via t\hat<tab>
ρ⃗ = f[i, CurvatureVector()]   # `ρ⃗` can be obtained via \rho<tab>\vec<tab>
@show t̂ ρ⃗
nothing  # hide

# We can check that both vectors are orthogonal and that their respective norms are 1 and
# approximately ``1/R`` (where ``R`` is the vortex ring radius):

using LinearAlgebra: norm, ⋅  # the dot product ̇`⋅` can be obtained via \cdot<tab>
@show t̂ ⋅ ρ⃗ norm(t̂) norm(ρ⃗)
nothing  # hide

# ### Plotting the filament

# We can readily plot our vortex ring using Makie.
# For convenience, VortexPasta overloads the Makie `plot` and `plot!` functions to be able
# to directly plot filaments.

using GLMakie
set_theme!(theme_black())

fig = Figure()                         # create an empty figure
ax = Axis3(fig[1, 1]; aspect = :data)  # add an Axis3 for plotting in 3D
zlims!(ax, 0.5, 1.5)                   # set axis limits in the z direction
plot!(ax, f)                           # plot filament onto axis
fig                                    # display the figure

# Note that, by default, the plot simply shows the `N` filament nodes (circular markers)
# joined by straight lines.
# To see how points are actually interpolated in-between nodes we can use the `refinement`
# keyword argument:

fig = Figure()
ax = Axis3(fig[1, 1]; aspect = :data)
zlims!(ax, 0.5, 1.5)
plot!(ax, f; refinement = 8)
fig

# It is also possible to plot other quantities such as the estimated tangent and curvature
# vectors along the curve:

fig = Figure()
ax = Axis3(fig[1, 1]; aspect = :data)
zlims!(ax, 0.5, 1.5)
plot!(
    ax, f;
    refinement = 8,
    tangents = true, tangentcolor = :Yellow,         # plot tangent vectors
    curvatures = true, curvaturecolor = :LightBlue,  # plot curvature vectors
    arrowscale = 0.16,  # rescale vectors
    vectorpos = 0.5,    # plot vectors at the midpoint in-between nodes
)
fig

# See [`filamentplot!`](@ref) for more details and for other possible options.

# ## Computing the vortex ring velocity
#
# An isolated vortex ring of radius ``R`` and circulation ``Γ`` is known to translate with a
# velocity [Saffman1993](@cite):
# ```math
# v_{\text{ring}} = \frac{Γ}{4πR} \left[ \ln \left(\frac{8R}{a}\right) - Δ \right],
# ```
# where ``a`` is the radius of the vortex core -- assumed to be much smaller than ``R`` --
# and ``Δ`` is a coefficient which depends on the actual vorticity profile at the core.
#
# !!! note "Some typical values of Δ"
#
#     - ``Δ = 1/2`` for a **hollow vortex**: ``ω(r) = ω₀ \, δ(r - a)``;
#     - ``Δ = 1/4`` for a **uniform vorticity distribution**: ``ω(r) = ω₀`` for ``r < a``.
#
# This velocity can be derived by computing the Biot--Savart integral along the circular
# vortex ring and excluding a very small region (proportional to ``a``) from the integral,
# in the vicinity of the point of interest, thus avoiding the singularity.
#
# In VortexPasta, the velocity induced by one or more vortex filaments is computed by the
# [`VortexPasta.BiotSavart`](@ref) submodule.
# The basic steps for computing the velocity induced by a set of vortices on itself is:
#
# 1. Set physical and numerical parameters ([`ParamsBiotSavart`](@ref)),
# 2. initialise a "cache" containing arrays and data needed for computations
#    ([`BiotSavart.init_cache`](@ref)),
# 3. compute filament velocities from their positions ([`velocity_on_nodes!`](@ref)).
#
# ### Physical and numerical parameters in periodic domains
#
# Relevant physical parameters are of two types:
#
# - vortex properties: circulation ``Γ``, core radius ``a`` and core parameter ``Δ``;
# - domain size ``L``, i.e. unit cell size in periodic domains.
#   Can be a single value (e.g. `Ls = 2π`) for a cubic domain, or a tuple of values (e.g.
#   `Ls = (2π, 4π, 6π)`) if one wants different domain sizes in each direction.
#
# These parameters are all mandatory.
# Moreover, there are a few numerical parameters as detailed in [`ParamsBiotSavart`](@ref).
# The two mandatory ones
