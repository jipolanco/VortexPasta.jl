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
# ```@contents
# Pages = ["01-vortex_ring.md"]
# Depth = 2:3
# ```
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
using VortexPasta.Filaments: Vec3

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

t̂ = f[i, UnitTangent()]       # "t̂" can be typed by t\hat<tab>
ρ⃗ = f[i, CurvatureVector()]   # "ρ⃗" can be typed by \rho<tab>\vec<tab>
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
    arrowwidth = 0.015, arrowscale = 0.14, arrowsize = (0.05, 0.05, 0.08),  # arrow properties
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
# ### Physical and numerical parameters
#
# Before being able to estimate the vortex ring velocity, we need to set the parameters
# needed to estimate Biot--Savart integrals.
# All of the required and optional parameters are described in [`ParamsBiotSavart`](@ref).
#
# Relevant **physical parameters** are of two types:
#
# - **vortex properties**: circulation ``Γ``, core radius ``a`` and core parameter ``Δ``;
# - **domain size** (or period) ``L``.

Γ = 1e-3   # vortex circulation                              || "Γ" can be typed by \Gamma<tab>
a = 1e-8   # vortex core size
Δ = 1/2    # vortex core parameter (1/2 for a hollow vortex) || "Δ" can be typed by \Delta<tab>
L = 2π     # domain period (same in all directions in this case)
nothing  # hide

# Note that here we work with nondimensional quantities, but we could try to relate these
# physical parameters to actual properties of, say, liquid helium-4.
# In this case, the vortex core size is ``a ≈ 10^{-8}\,\text{cm}`` and the quantum of
# circulation is ``Γ = κ = h/m ≈ 0.997 × 10^{-3}\,\text{cm}^2/\text{s}``, where ``h`` is
# Planck's constant and ``m`` the mass of a helium atom.
# So the above parameters are quite relevant to ``^4\text{He}`` if one interprets lengths
# in centimetres and times in seconds.
#
# There are also important **numerical parameters** which need to be set.
# See the [Parameter selection](@ref Ewald-parameters) section for an advice on how to set
# them.
#
# Following that page, we start by setting the resolution ``N`` of the numerical grid used
# for long-range computations, and we set the other parameters based on that:

N = round(Int, 64 * 4/5)  # we prefer if the FFT size is a power of 2, here Ñ = σN = 64
kmax = π * N / L          # this is the maximum resolved wavenumber (the Nyquist frequency)
α = kmax / 8              # Ewald splitting parameter || "α" can be typed by "\alpha<tab>"
rcut = 5 / α              # cut-off distance for short-range computations
rcut / L                  # note: the cut-off distance should be less than half the period L

# Additionally, we can optionally set the parameters for numerical integration.
# These are the quadrature rules used to approximate line integrals within each filament
# segment (see [Numerical integration](@ref)).
# We can separately set the quadrature rule used for short-range and long-range computations:

quadrature_short = GaussLegendre(4)  # use 4-point Gauss–Legendre quadrature
quadrature_long  = GaussLegendre(4)
nothing  # hide

# See [`GaussLegendre`](@ref) or the [Wikipedia
# page](https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss%E2%80%93Legendre_quadrature) for more details.
#
# Finally, we put all these parameters together in a [`ParamsBiotSavart`](@ref) object:

using VortexPasta.BiotSavart
params = ParamsBiotSavart(;
    Γ, α, Δ, a,
    Ls = (L, L, L),  # same domain size in all directions
    Ns = (N, N, N),  # same long-range resolution in all directions
    rcut,
    quadrature_short,
    quadrature_long,
)

# Note that there are a few parameters, namely the short-range and long-range backends,
# which haven't been discussed yet.
# These correspond to different ways of computing both components.
# We leave them here at their default values.
# See [`CellListsBackend`](@ref) and [`FINUFFTBackend`](@ref) for details on the defaults
# and their parameters.

# ### Computing the velocity on filament nodes
#
# The last thing to do before launching computations is to create a "cache" where all data
# (arrays, …) required for computations are kept.
# This is in part to avoid making large allocations every time we wish to compute
# Biot--Savart integrals, which would be quite bad for performance.
# Luckily, creating a cache is very simple:

fs = [f]  # note: we need to work with a *vector* of filaments
cache = BiotSavart.init_cache(params, fs)
nothing  # hide

# We can now estimate the velocity of all filament nodes using [`velocity_on_nodes!`](@ref).
# But first we need to allocate the output, which will contain all velocities:

vs = map(similar ∘ nodes, fs)
@show summary(vs) summary(vs[1])
nothing  # hide

# Note that this is a 1-element vector (because we only have 1 filament).
# Moreover, the element `vs[1]` is a vector of [`Vec3`](@ref) (3-element vectors,
# specifically velocity vectors in this case), in which the element
# `vs[1][i]` will be the velocity of the filament node `fs[1][i]`.
#
# To actually fill `vs` with the velocities of the filament nodes, we need to call
# [`velocity_on_nodes!`](@ref):

velocity_on_nodes!(vs, cache, fs)
vs[1]  # prints the velocities of all nodes of the first (and only) filament

# We have done our first Biot--Savart calculation!
# We can see that all nodes have practically zero velocity in the ``x`` and ``y``
# directions, while they all have approximately the same velocity in the ``z`` direction.
# Good news, this is exactly what we expect for a vortex ring!
#
# We can compute the mean velocity and quantify the velocity variation along the filament:

using Statistics: mean, std
v_mean = mean(vs[1])
v_std = std(vs[1])
v_std_normalised = v_std / norm(v_mean)  # normalised standard deviation
@show v_mean v_std v_std_normalised
nothing  # hide

# Note that there are some (very small) fluctuations of the velocity about the mean.
# As discussed and verified further below, this is explained by periodicity effects, as the
# vortex ring is affected by its periodic images.
# Furthermore, different points along the ring are at different distances from the images,
# and therefore the velocity induced by the images varies slightly as we move along the
# vortex ring.
# In other words, periodicity introduces a slight anisotropy on the behaviour of the vortex ring.
#
# To visualise things, it is easy to plot the filament with the node velocities:

fig = Figure()
ax = Axis3(fig[1, 1]; aspect = :data)
zlims!(ax, 0.5, 1.5)
plot!(
    ax, fs[1], vs[1];
    refinement = 4, arrowsize = (0.03, 0.03, 0.04), arrowwidth = 0.01, arrowscale = 30,
)
fig

# ### Comparison with analytical expression
#
# We can finally compare the filament velocity with the analytical vortex ring velocity
# according to the expression at the [start of this section](@ref
# Computing-the-vortex-ring-velocity):

v_ring = Γ / (4π * R) * (log(8R / a) - Δ)
vz = v_mean[3]
relative_difference = (v_ring - vz) / v_ring
@show v_ring vz relative_difference
nothing  # hide

# We can see that we're very close, and that the relative difference between the two is of
# less than 0.1%.
# This is already very nice, but note that some of the difference may be explained by
# **periodicity effects**.
# Indeed, the vortex ring is not alone, but it is also affected by its periodic images.
#
# We can visualise this by plotting the ring and its nearest 26 periodic images (actually
# there's an infinity of them!):

fig = Figure()
ax = Axis3(fig[1, 1]; aspect = :data)
hidespines!(ax)
hidexdecorations!(ax; label = false)
hideydecorations!(ax; label = false)
hidezdecorations!(ax; label = false)
fimage = copy(f)
plot!(ax, f; refinement = 4, markersize = 0, color = :OrangeRed, linewidth = 2)
for I ∈ CartesianIndices((-1:1, -1:1, -1:1))
    widths = (L, L, L)
    offset = Vec3(Tuple(I) .* widths)
    box = Rect(offset..., widths...)
    wireframe!(ax, box; color = (:white, 0.5), linewidth = 0.2)  # plot cube
    iszero(I) && continue  # don't replot the "original" filament
    fimage .= f .+ Ref(offset)
    update_coefficients!(fimage)  # updates interpolation coefficients
    plot!(ax, fimage; refinement = 4, markersize = 0, color = :LightBlue, linewidth = 1)
end
fig

# One can show that the images have the tendency to slow the vortex down, which is
# consistent with our results (since `vz < v_ring`).
# Of course, the effect here is negligible, but it should become more noticeable for larger
# vortex rings (i.e. with diameter ``2R`` close to the domain size ``L``).

# ### Side note: disabling periodicity
#
# It is actually possible to disable the periodic boundary conditions, effectively leading
# to an open domain with a unique vortex ring.
#
# In VortexPasta.jl, disabling periodicity amounts to disabling the long-range part of Ewald
# summation, so that **all** scales are computed by the "short"-range component.
# The cut-off distance is therefore infinity, and it makes no sense to use cell lists to
# accelerate the computations.
# Other methods exist (and are now standard) to accelerate this kind of computation, but
# these are not implemented here.
# So, while it is possible to disable periodicity, this should only be used for testing or
# for small calculations.
#
# To disable periodicity, one should pass `α = Zero()` and `Ls = Infinity()` to
# [`ParamsBiotSavart`](@ref):

params_inf = ParamsBiotSavart(;
    Γ, α = Zero(), Δ, a,
    Ls = Infinity(),
    quadrature_short,
    quadrature_long,
)

# Here [`Zero`](@ref) and [`Infinity`](@ref) are custom types that, as one may expect,
# represent ``0`` and ``+∞``.
#
# We can now compute the vortex ring velocity in the absence of periodic effects:

cache_inf = BiotSavart.init_cache(params_inf, fs)
vs_inf = map(similar, vs)
velocity_on_nodes!(vs_inf, cache_inf, fs)
v_mean_inf = mean(vs_inf[1])
v_std_inf = std(vs_inf[1])
v_std_normalised_inf = v_std_inf / norm(v_mean_inf)
@show v_mean_inf v_std_inf v_std_normalised_inf
nothing  # hide

# Interestingly, the velocity fluctuations along the filament are orders of magnitude
# smaller than in the periodic case.
# This is consistent with what was explained on the anisotropy introduced by periodicity,
# in the sense that not all points on the filament "feel" the effect of the periodic images
# in the same way.
#
# Besides, the mean velocity is approximately the same as in the periodic case.
# In fact, we can see that it's even closer to the analytical velocity:

vz_inf = v_mean_inf[3]
relative_difference_inf = (v_ring - vz_inf) / v_ring

# The relative difference is close to 0.01%!
# Once again, this clearly shows the (negligible) effect of periodic images on the effective
# vortex ring velocity.

# ## Making the ring move
#
# Up to now we have shown how to compute the velocity self-induced by a vortex filament on
# its discretisation points.
# We now want to use this velocity to advect the vortex filament over time.
# In the case of a vortex ring, its displacement is very simple and boring (the vortex ring
# moves with constant velocity), but this case is still useful for instance for testing the
# stability of timestepping schemes.
#
# For running temporal simulations, we want to use the [`VortexPasta.Timestepping`](@ref)
# submodule.
# The syntax is somewhat inspired from the popular
# [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) ecosystem.
#
# We start by defining a [`VortexFilamentProblem`](@ref), which takes (1) an initial condition
# (in our case, the initial vortex ring we already defined); (2) the time span of the
# simulation; and (3) parameters for Biot--Savart computations (which we already defined as well):

using VortexPasta.Timestepping
T = L / 2vz       # time it should take for the ring to cross half the periodic box
tspan = (0.0, T)  # time span; this basically determines when to stop the simulation
prob = VortexFilamentProblem(fs, tspan, params)

# The second step is to "initialise" the problem, which means choosing a [temporal
# scheme](@ref Temporal-schemes), an initial timestep, as well as optional parameters
# related for instance to [temporal adaptivity](@ref Adaptivity), [spatial refinement](@ref
# Refinement) of the filaments, or [vortex reconnections](@ref Reconnections).
#
# For now, we use the explicit [`RK4`](@ref) timestepping scheme, set a constant timestep
# `dt`, and leave the defaults for everything else:

l_min = minimum_knot_increment(fs)   # this is an estimate for the minimum distance between discretisation points
dt = 20 * l_min^2 / (Γ * log(l_min / a))
iter = init(prob, RK4(); dt)

# Note that the maximum possible timestep to avoid numerical instability is mainly limited
# by the filament resolution ``ℓ_{\text{min}}``.
# More precisely, the limit is the frequency of spurious oscillations associated to that
# length scale.
# These correspond to Kelvin waves, whose frequency roughly scales as
# ``ω(ℓ) ∼ Γ \ln (ℓ / a) / ℓ^2``.
#
# We now run the simulation until ``t = T = L / 2 v_{\text{ring}}``:

solve!(iter)      # run the simulation until t = T
iter.time.t / T  # check that the current time is t == T

# We can check that, as expected, the ring has crossed half the periodic box:

displacement = @. (iter.fs - prob.fs) / L
displacement[1]  # prints displacement of the first (and only) filament

# As we can see, the filament nodes have basically only moved in the ``z`` direction by
# almost exactly ``L / 2``.
