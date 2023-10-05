# # [Kelvin waves](@id tutorial-kelvin-waves)
#
# This tutorial describes the simulation of Kelvin waves propagating along nearly-straight
# and infinite vortex lines.
#
# Here we will:
#
# - learn how to define infinite but unclosed filaments;
# - look at diagnostics such as the energy over time.
#
# It is recommended to first follow the [vortex ring tutorial](@ref tutorial-vortex-ring)
# before following this tutorial.
#
# ```@contents
# Pages = ["02-kelvin_waves.md"]
# Depth = 2:3
# ```
#
# ## Initialising infinite lines
#
# The idea of this tutorial is to study the time evolution of an infinite straight line
# slightly modified by a sinusoidal perturbation.
#
# We will consider such a vortex line in a cubic periodic domain of size ``L = 2π``.
# The line is oriented in the ``z`` direction and modified by a perturbation of
# amplitude ``ϵL`` along ``x``.
# The perturbation is periodic with period ``λ = L/m = 2π/m`` where ``m`` is
# a positive integer.
#
# Such an infinite line ``\bm{s} = (x, y, z)`` can be parametrised as
#
# ```math
# \begin{align*}
#   x(t) &= x_0 + ϵ \, \sin(2πmt) \\
#   y(t) &= y_0 \\
#   z(t) &= z_0 + 2π (t - 1/2)
# \end{align*}
# ```
#
# for ``t ∈ \mathbb{R}``.
# In particular, note that the line exactly jumps from ``\bm{s}(t)`` to
# ``\bm{s}(t + 1) = \bm{s}(t) + (0, \, 0, \, 2π)`` after a period ``T = 1``.
#
# ### Defining an unclosed infinite curve
#
# Following the [vortex ring tutorial](@ref tutorial-vortex-ring-init-filament), one may
# want to define such a line as follows:

using VortexPasta
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3

N = 16  # number of discretisation points per line
m = 2   # perturbation mode
L = 2π
x⃗₀ = Vec3(L/4, L/4, L/2)
ϵ = 0.02
ts = range(0, 1; length = N + 1)[1:N]
points = [x⃗₀ + Vec3(ϵ * L * sinpi(2m * t), 0, L * (t - 1/2)) for t ∈ ts]
f = Filaments.init(ClosedFilament, points, CubicSplineMethod())
nothing  # hide

# Let's look at the result:

using GLMakie
set_theme!()  # hide
GLMakie.activate!()  # hide

## Plot a list of filaments
function plot_filaments(fs::AbstractVector{<:AbstractFilament})
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    hidespines!(ax)
    wireframe!(ax, Rect(0, 0, 0, L, L, L); color = (:black, 0.5), linewidth = 0.2)
    for f ∈ fs
        plot!(ax, f; refinement = 4)
    end
    fig
end

## Plot a single filament
plot_filaments(f::AbstractFilament) = plot_filaments([f])

plot_filaments(f)

# Things look almost as expected except for the fact that the line tries to close itself
# when it reaches the end.
# To avoid this, one needs to explicitly give [`Filaments.init`](@ref) an end-to-end
# vector via the `offset` keyword argument.
# In our case the end-to-end vector is ``\bm{Δ} = (0, 0, 2π)``.

f = Filaments.init(ClosedFilament, points, CubicSplineMethod(); offset = (0, 0, 2π))
plot_filaments(f)

# Now everything looks fine!
# Note that the end-to-end vector corresponds to the separation between a node `f[i]` and
# the node `f[i + N]`.
# For example:

@show f[end + 1] - f[begin]
@show Vec3(0, 0, 2π)
nothing  # hide

# !!! note "End-to-end vector"
#     
#     The end-to-end vector *must* be an integer multiple of the domain period, which in
#     this case is ``(2π, 2π, 2π)``.

# ### Using predefined curves
#
# There is an alternative way of defining such curves, using the
# [`VortexPasta.PredefinedCurves`](@ref) module which provides convenient definitions of
# commonly-used curves.
#
# In this case we want to use the [`PeriodicLine`](@ref) definitions, which allow one to
# pass arbitrary functions as perturbations.
# Note that curve definitions in `PredefinedCurves` are normalised.
# In particular, the period of `PeriodicLine` is 1, and the perturbation that we give it
# will be in terms of this unit period.

using VortexPasta.PredefinedCurves: PeriodicLine, define_curve
x_perturb(t) = ϵ * sinpi(2m * t)  # perturbation along x (takes t ∈ [0, 1])
p = PeriodicLine(x = x_perturb)      # this represents a line with period 1 along z
nothing  # hide

# We now want to "convert" this line to a parametric function which can be then evaluated to
# generate points.
# This is done using the [`define_curve`](@ref) function, which allows in particular to
# rescale the curve (we want a period of ``L = 2π`` instead of ``1``).
# We would also like to shift the origin to ``\bm{x}_0``.

S = define_curve(p; scale = L, translate = x⃗₀)
@show S(0.0) S(0.5) S(1.0)
nothing  # hide

# As we can see, `S` is a function which can be evaluated at any value of ``t``.
# We could generate filament nodes by evaluating it over a range of ``t``s, or we can let
# [`Filaments.init`](@ref) do that for us:

f′ = Filaments.init(S, ClosedFilament, N, CubicSplineMethod())
nothing  # hide

# Note that this generates a filament which is practically identical to the previous one
# (just with a shift in the node positions):

plot_filaments([f, f′])

# ### Ensuring periodicity of the velocity
#
# For now we have initialised one infinite unclosed filament.
# One needs to be careful when working with unclosed filaments in periodic domains.
# Indeed, a single straight vortex filament in a periodic domain generates a non-zero
# circulation along the domain boundaries (or equivalently, a non-zero mean vorticity),
# which violates the periodicity condition.
#
# Note that the mean vorticity in the periodic domain is given by
#
# ```math
# ⟨ \bm{ω} ⟩
# = \frac{1}{V} ∫_Ω \bm{ω}(\bm{x}) \, \mathrm{d}^3\bm{x}
# = \frac{Γ}{V} ∫_{\mathcal{C}} \mathrm{d}\bm{s}
# = \frac{Γ}{V} ∫_{\mathcal{C}} \bm{s}'(t) \, \mathrm{d}t
# ```
#
# where ``Ω`` represents the periodic domain and ``V = L^3`` is its volume.
# So the last integral must be zero to ensure periodicity.
# It is quite obvious that this is not the case for the filament defined above,
# and we can readily verify it:

integrate(f′, GaussLegendre(4)) do ff, i, ζ
    ff(i, ζ, Derivative(1))
end

# This means that, for each vortex oriented in the ``+z`` direction, we need to compensate
# by a vortex oriented in the ``-z`` direction to obtain a zero total circulation.
#
# Secondly, if we put just two vortices of opposite signs in the domain, these will rotate
# and "dance" around each other due to the velocity they induce on each other.
# To avoid this and stabilise the system, one can add two more vortices, in such a way that
# all vortices are equally spaced on the XY plane.
# Let's create these four vortices:

funcs = [
    # "Positive" vortices
    define_curve(p; scale = L, translate = (0.25L, 0.25L, 0.5L)),
    define_curve(p; scale = L, translate = (0.75L, 0.75L, 0.5L)),
    # "Negative" vortices: we scale by -L to change the line orientation
    define_curve(p; scale = -L, translate = (0.25L, 0.75L, 0.5L)),
    define_curve(p; scale = -L, translate = (0.75L, 0.25L, 0.5L)),
]
fs = map(S -> Filaments.init(S, ClosedFilament, N, CubicSplineMethod()), funcs)
plot_filaments(fs)

# We can check that, when we sum the contributions of all filaments, the mean vorticity is
# zero:

sum(fs) do f
    integrate(f, GaussLegendre(4)) do ff, i, ζ
        f(i, ζ, Derivative(1))
    end
end

# Now we're ready to perform simulations.
