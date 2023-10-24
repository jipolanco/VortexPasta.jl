# # [Kelvin waves](@id tutorial-kelvin-waves)
#
# This tutorial describes the simulation of Kelvin waves propagating along nearly-straight
# and infinite vortex lines.
#
# Here we will:
#
# - learn how to define infinite but unclosed filaments;
# - look at diagnostics such as the energy over time;
# - perform spatial and temporal Fourier analysis to detect relevant wavenumbers and
#   frequencies associated to Kelvin waves.
#
# It is recommended to first follow the [vortex ring tutorial](@ref tutorial-vortex-ring)
# before following this tutorial.
#
# ```@contents
# Pages = ["02-kelvin_waves.md"]
# Depth = 2:3
# ```
#
# ## Physical configuration
#
# The idea of this tutorial is to study the time evolution of an infinite straight line
# slightly modified by a sinusoidal perturbation.
#
# We will consider such a vortex line in a cubic periodic domain of size ``L = 2π``.
# The line is oriented in the ``z`` direction and modified by a perturbation of
# amplitude ``ϵL`` along ``x``.
# The perturbation is periodic with period ``λ = L/m = 2π/m`` where ``m`` is an integer
# representing the *mode* of the perturbation (relative to the domain size ``L``).
#
# Such an infinite line ``\bm{s} = (x, y, z)`` can be parametrised as
#
# ```math
# \begin{align*}
#   x(t) &= x_0 + ϵ \, \sin(2πmt) \\
#   y(t) &= y_0 \\
#   z(t) &= z_0 + \left( t - \frac{1}{2} \right) L
# \end{align*}
# ```
#
# for ``t ∈ \mathbb{R}``.
# In particular, note that the line exactly crosses the domain after a period ``T = 1``,
# going from ``\bm{s}(t)`` to ``\bm{s}(t + 1) = \bm{s}(t) + (0, \, 0, \, L)``.
#
# The analytical prediction is that, over time, a small perturbation should rotate around
# the vortex in the direction opposite to its circulation.
# Its frequency is given by (see e.g. [Schwarz1985](@citet)):
#
# ```math
# ω_{\text{KW}}(k) = \frac{Γ k^2}{4π} \left[
#   \ln\left( \frac{2}{k a} \right) - γ + \frac{1}{2} - Δ
# \right]
# ```
#
# where ``k = 2πm/L`` is the perturbation wavenumber, ``γ ≈ 0.5772`` the [Euler–Mascheroni
# constant](https://en.wikipedia.org/wiki/Euler%27s_constant), ``Δ`` the vortex core
# parameter and ``a`` its radius (``a ≪ 1/k``).
#
# ## Defining an unclosed infinite curve
#
# Following the [vortex ring tutorial](@ref tutorial-vortex-ring-init-filament), one may
# want to define such a line as follows:

using VortexPasta
using VortexPasta.Filaments
using VortexPasta.Filaments: Vec3

N = 64  # number of discretisation points per line
m = 2    # perturbation mode
L = 2π   # domain period
x⃗₀ = Vec3(L/4, L/4, L/2)  # line "origin"
ϵ = 0.01
ts = range(0, 1; length = N + 1)[1:N]  # important: we exclude the endpoint (t = 1)
points = [x⃗₀ + Vec3(ϵ * L * sinpi(2m * t), 0, L * (t - 1/2)) for t ∈ ts]
f = Filaments.init(ClosedFilament, points, CubicSplineMethod())
nothing  # hide

# Let's look at the result:

using GLMakie
set_theme!()  # hide
GLMakie.activate!()  # hide

## Give a colour to a filament based on its local orientation wrt Z.
function filament_colour(f::AbstractFilament, refinement)
    cs = Float32[]
    ζs = range(0, 1; length = refinement + 1)[1:refinement]  # for interpolation
    for seg ∈ segments(f), ζ ∈ ζs
        colour = seg(ζ, UnitTangent())[3]  # in [-1, 1]
        push!(cs, colour)
    end
    cs
end

## Plot a list of filaments
function plot_filaments(fs::AbstractVector)
    fig = Figure()
    ax = Axis3(fig[1, 1]; aspect = :data)
    ticks = range(0, 2π; step = π/2)
    tickformat(xs) = map(x -> string(x/π, "π"), xs)
    ax.xticks = ax.yticks = ax.zticks = ticks
    ax.xtickformat = ax.ytickformat = ax.ztickformat = tickformat
    hidespines!(ax)
    wireframe!(ax, Rect(0, 0, 0, L, L, L); color = (:black, 0.5), linewidth = 0.2)
    for f ∈ fs
        refinement = 4
        color = filament_colour(f, refinement)
        plot!(
            ax, f;
            refinement, color, colormap = :RdBu_9, colorrange = (-1, 1), markersize = 4,
        )
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

# ## Defining a curve from parametric function
#
# The [`Filaments.init`](@ref) function actually allows to define a curve directly from its
# (continuous) parametric function.
# In this case one doesn't need to care about end-to-end vectors and "offsets", since these
# are usually encoded in the parametrisation.
#
# For example, for the curve above we would define the function:

fcurve(t) = x⃗₀ + Vec3(
    ϵ * L * sinpi(2 * m * t),
    0,
    (t - 0.5) * L,
)

# The function will be evaluated over the interval ``t ∈ [0, 1]``.
# The only assumption is that the parametric function must either represent:
#
# - a closed curve with period ``T = 1``;
# - an unclosed periodic curve which crosses the domain after a period ``T = 1``.
#
# Here we are in the second case, and the function above indeed satisfies this condition.
#
# Now we just pass the function to `Filaments.init`:

f′ = Filaments.init(fcurve, ClosedFilament, N, CubicSplineMethod())
nothing  # hide

# Note that this generates a filament which is practically identical to the previous one
# (just with a shift in the node positions, not really visible here):

plot_filaments([f, f′])

# ## Using predefined curves
#
# There is another convenient way of defining such curves, using the
# [`VortexPasta.PredefinedCurves`](@ref) module which provides definitions of parametric
# functions for commonly-used curves.
# As we will see in the next section, this is particularly convenient when we want to
# create multiple vortices which share the same geometry, but which have for instance
# different orientations or different spatial locations in the domain.
#
# Here we want to use the [`PeriodicLine`](@ref) definitions, which allow one to
# pass arbitrary functions as perturbations.
# Note that curve definitions in `PredefinedCurves` are normalised.
# In particular, the period of `PeriodicLine` is 1, and the perturbation that we give it
# will be in terms of this unit period.

using VortexPasta.PredefinedCurves: PeriodicLine, define_curve
x_perturb(t) = ϵ * sinpi(2m * t)  # perturbation along x (takes t ∈ [0, 1])
p = PeriodicLine(x = x_perturb)   # this represents a line with period 1 along z
nothing  # hide

# We now want to "convert" this line to a parametric function which can be then evaluated to
# generate points.
# This is done using the [`define_curve`](@ref) function, which allows in particular to
# rescale the curve (we want a period of ``L = 2π`` instead of ``1``).
# We would also like the curve to be centred at ``\bm{x}_0``.

S = define_curve(p; scale = L, translate = x⃗₀)
@show S(0.0) S(0.5) S(1.0)
nothing  # hide

# As we can see, `S` is a function which can be evaluated at any value of ``t``.
# In fact, `S` is identical to the `fcurve` function we defined above.
# We can now pass this function to [`Filaments.init`](@ref) to generate a filament:

f = Filaments.init(S, ClosedFilament, N, CubicSplineMethod())
plot_filaments(f)

# ## Ensuring periodicity of the velocity
#
# For now we have initialised one infinite unclosed filament.
# One needs to be careful when working with unclosed filaments in periodic domains.
# Indeed, a single straight vortex filament in a periodic domain generates a non-zero
# circulation along the domain boundaries (or equivalently, a non-zero mean vorticity),
# which violates the periodicity condition.
#
# The mean vorticity in the periodic domain is given by
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
# To avoid this and stabilise the system, we can add two more vortices, making sure that
# all vortices are equally spaced on the XY plane.
# Let's create these four vortices:

funcs = [
    ## "Positive" vortices
    define_curve(p; scale = L, translate = (0.25L, 0.25L, 0.5L)),
    define_curve(p; scale = L, translate = (0.75L, 0.75L, 0.5L)),
    ## "Negative" vortices: we use the `orientation` keyword to flip their orientation.
    define_curve(p; scale = L, translate = (0.25L, 0.75L, 0.5L), orientation = -1),
    define_curve(p; scale = L, translate = (0.75L, 0.25L, 0.5L), orientation = -1),
]
fs = map(S -> Filaments.init(S, ClosedFilament, N, CubicSplineMethod()), funcs)
plot_filaments(fs)

# Here the colours represent the local orientation of the curve tangent with respect to the
# ``z`` axis.
# We can check that, when we sum the contributions of all filaments, the mean vorticity is
# zero:

## This computes the integral along each filament and sums the results.
sum(fs) do f
    integrate(f, GaussLegendre(4)) do ff, i, ζ
        f(i, ζ, Derivative(1))
    end
end

# Now we're ready to perform simulations.

# ## Simulating Kelvin waves
#
# As in the [vortex ring tutorial](@ref tutorial-vortex-ring), we use the
# [`Timestepping`](@ref) module to perform a temporal simulation of the configuration we
# just prepared.
#
# ### Setting physical and numerical parameters
#
# We start by setting the parameters for Biot--Savart computations:

using VortexPasta.BiotSavart
M = round(Int, 32 * 4/5)  # resolution of long-range grid
kmax = π * M / L          # maximum resolved wavenumber (Nyquist frequency) for long-range part
α = kmax / 8              # Ewald splitting parameter

params = ParamsBiotSavart(;
    Γ = 1.0,    # vortex circulation
    a = 1e-8,   # vortex core size
    Δ = 1/4,    # vortex core parameter (1/4 for a constant vorticity distribution)
    α = α,      # Ewald splitting parameter
    Ls = (L, L, L),  # same domain size in all directions
    Ns = (M, M, M),  # same long-range resolution in all directions
    rcut = 5 / α,    # cut-off distance for short-range computations
    quadrature_short = GaussLegendre(2),  # quadrature for short-range computations
    quadrature_long = GaussLegendre(2),   # quadrature for long-range computations
    backend_long = FINUFFTBackend(),      # this is the default
    backend_short = NaiveShortRangeBackend(),  # OK when rcut/L is large
)

# We would like to compute a few periods of Kelvin wave oscillations.
# For this, we first compute the expected Kelvin wave frequency and its associated period:

(; Γ, a, Δ,) = params         # extract parameters needed for KW frequency
γ = MathConstants.eulergamma  # Euler–Mascheroni constant
k = 2π * m / L
ω_kw = Γ * k^2 / (4 * π) * (
    log(2 / (k * a)) - γ + 1/2 - Δ
)
T_kw = 2π / ω_kw              # expected Kelvin wave period

# We create a [`VortexFilamentProblem`](@ref) to simulate a few Kelvin wave periods:

using VortexPasta.Timestepping
tspan = (0.0, 4 * T_kw)
prob = VortexFilamentProblem(fs, tspan, params)

# We now create a callback which will be used to store some data for further analysis.
# We will store the times and the position over time of a single filament node to be able to
# visualise and analyse the oscillations.
# Moreover, we will store the system energy to verify that energy is conserved over time
# (see [VFM notes](@ref VFM-energy) for detains on how it is computed).
# For computing the energy we use the [`kinetic_energy_from_streamfunction`](@ref) function
# from the [`Diagnostics`](@ref) module.

using VortexPasta.Diagnostics

times::Vector{Float64} = Float64[]
X_probe::Vector{Vec3{Float64}} = Vec3{Float64}[]  # will contain positions of a chosen node
energy::Vector{Float64} = Float64[]

function callback(iter)
    (; nstep, t,) = iter
    if nstep == 0  # make sure vectors are empty at the beginning of the simulation
        empty!(times)
        empty!(X_probe)
        empty!(energy)
    end
    push!(times, t)
    s⃗ = iter.fs[1][3]  # we choose a single node of a single filament
    push!(X_probe, s⃗)
    ## Compute energy
    E = Diagnostics.kinetic_energy_from_streamfunction(iter)
    push!(energy, E)
    nothing
end

# Note that we have annotated the types of the variables `times`, `X_probe` and `energy` for
# performance reasons, since these are *global* variables which are used (and modified) from
# within the `callback` function.
# See
# [here](https://docs.julialang.org/en/v1/manual/performance-tips/#Avoid-untyped-global-variables)
# and
# [here](https://docs.julialang.org/en/v1/manual/variables-and-scoping/#man-typed-globals)
# for details.

# ### Choosing the timestep and the temporal scheme
#
# In the [vortex ring tutorial](@ref tutorial-vortex-ring) we have used the standard
# [`RK4`](@ref) scheme.
# To capture the vortex evolution and avoid blow-up, this scheme requires the timestep
# ``Δt`` to be of the order of the period of the **fastest** resolved Kelvin waves, which have a
# wavelength ``λ`` equal to (twice) the filament resolution ``δ`` (the typical distance
# between two discretisation points).
#
# We first estimate the filament resolution using [`minimun_node_distance`](@ref):

δ = minimum_node_distance(prob.fs)  # should be close to L/N in our case

# Now we compute the Kelvin wave frequency associated to this distance:

kelvin_wave_period(λ; a, Δ, Γ) = 2 * λ^2 / Γ / (log(λ / (π * a)) + 1/2 - (Δ + MathConstants.γ))
dt_kw = kelvin_wave_period(δ; a, Δ, Γ)

# Note that this time scale is very small compared to the period of the large-scale Kelvin
# waves:

T_kw / dt_kw

# This means that we would need a relatively large simulation time to observe the
# evolution of large-scale Kelvin waves over multiple periods.

# Besides, for the [`RK4`](@ref) scheme, this time scale really seems to set the maximum allowed
# timestep limit.
# We can check that a simulation with `RK4` using `dt = dt_kw` remains stable.
# In particular, energy stays constant over time after running a few iterations:

iter = init(prob, RK4(); dt = dt_kw, callback)
for _ ∈ 1:40
    step!(iter)
end
energy'

# Meanwhile, using `dt = 2 * dt_kw` quickly leads to instability and energy blow-up:

iter = init(prob, RK4(); dt = 2 * dt_kw, callback)
for _ ∈ 1:40
    step!(iter)
end
energy'

# This limit is basically set by the local self-induced velocity term (see
# [Desingularisation](@ref VFM-desingularisation)), which presents fast temporal variations.
# This term is actually very cheap to compute compared to the non-local velocity resulting
# from the Biot--Savart law.
# In this case, it makes sense to use multirate timestepping schemes, which use a small
# timestep to evaluate the fast (and cheap) motions, and a larger timestep to evaluate the
# slow (and expensive) terms.
# In the following we use the [`SanduMRI33a`](@ref) scheme [Sandu2019](@cite), which is a
# 3rd-order 3-stage Runge--Kutta scheme for the "slow" component.
# For each "outer" RK stage, we perform ``M = 10`` inner iterations of the "fast" component
# using the standard `RK4` scheme:

scheme = SanduMRI33a(RK4(), 10)

# This allows to greatly increase the timestep:

dt = 32 * dt_kw  # timestep allowed by SanduMRI33a(RK4(), 10) scheme

# Note that we could tune the number ``M`` of inner iterations to allow even larger
# timesteps, but this can lead to precision loss (in particular, small energy fluctuations).
#
# ### Running the simulation
#
# We now solve the full problem with this multirate scheme:

iter = init(prob, scheme; dt, callback)
reset_timer!(iter.to)  # to get more accurate timings (removes most of the compilation time)
solve!(iter)

# We can check that energy is conserved:

energy'

# We see that the energy seems to take the same value at all times.
# We can verify this quantitatively by looking at its standard deviation (normalised by the
# mean energy), which is negligible:

using Statistics: mean, std
Emean = mean(energy)
Estd = std(energy)
Estd / Emean

# We now plot the evolution of the ``x`` and ``y`` coordinates of the closen filament node:

using CairoMakie  # hide
CairoMakie.activate!(type = "svg", pt_per_unit = 1.0)  # hide
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"t / T_{\text{KW}}", ylabel = "Position")
tnorm = times ./ T_kw  # normalised time
xpos = map(s⃗ -> s⃗[1], X_probe)  # get all X positions over time
ypos = map(s⃗ -> s⃗[2], X_probe)  # get all Y positions over time
scatterlines!(ax, tnorm, xpos; marker = 'o', label = L"x(t)")
scatterlines!(ax, tnorm, ypos; marker = 'o', label = L"y(t)")
Legend(fig[1, 2], ax; orientation = :vertical, framevisible = false, padding = (0, 0, 0, 0), labelsize = 22, rowgap = 5)
fig

# We see that the ``x`` and ``y`` positions of the chosen point oscillate sinusoidally.
# The period of the oscillations are very close to the expected Kelvin wave period
# ``T_{\text{KW}}``.
#
# The oscillations above suggest circular trajectories, as we can check in the following
# figure:

scatterlines(
    xpos, ypos;
    color = tnorm,
    axis = (aspect = DataAspect(), xlabel = L"x(t)", ylabel = L"y(t)"),
)

# ## Measuring performance
#
# The VortexPasta solver uses the
# [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl) package
# to estimate the time spent (and memory allocated) in different stages of the computation.
#
# Accessing timing information is very simple, as it is all included in the `to` field of
# the [`VortexFilamentSolver`](@ref Timestepping.VortexFilamentSolver):

show(IOContext(stdout, :displaysize => (40, 100)), iter.to)  # hide
iter.to
nothing  # hide

# We can see that, in this case, roughly half the time is spent in the long-range
# computations, while the other half is spent on short-range computations and the LIA
# (local) term.
# Note that the LIA term is computed much more often than the other components due to the
# use of a multirate scheme, but even then, its total cost is a fraction of that associated
# to the Biot--Savart integrals.

# ## Fourier analysis
#
# ### Spatial analysis
#
# The idea is to identify the spatial fluctuations of a single vortex with respect to the
# unperturbed filament.
# For this, we first write the perturbations in complex representation as a function of the ``z``
# coordinate, i.e. ``r(z) = x(z) + i y(z)``.
#
# We want to apply the FFT to these two functions.
# For this, we need all points of the vortex filament to be equispaced in ``z``:

f = iter.fs[1]               # vortex to analyse
zs = getindex.(nodes(f), 3)  # z locations
N = length(zs)
zs_expected = range(zs[begin], zs[begin] + L; length = N + 1)[1:N]  # equispaced locations
isapprox(zs, zs_expected; rtol = 1e-5)  # check that z locations are approximately equispaced

# Now that we have verified this, we define the complex function ``r(z) = x(z) + i y(z)``
# and we perform a complex-to-complex FFT to obtain ``\hat{r}(k)``:

xs = getindex.(nodes(f), 1)  # x locations
ys = getindex.(nodes(f), 2)  # y locations
rs = @. xs + im * ys

using FFTW: fft, fft!, fftfreq
rhat = fft(rs)
@. rhat = rhat / N  # normalise FFT
@show rhat[1]       # the zero frequency gives the mean location
rhat[1] ≈ π/2 + π/2 * im  # we expect the mean location to be (π/2, π/2)

# The associated wavenumbers are multiples of ``2π/Δz = 2πN/L``:

Δz = L / N
@assert isapprox(Δz, zs[2] - zs[1]; rtol = 1e-4)
ks = fftfreq(N, 2π / Δz)
ks'  # should be integers if L = 2π

# Note that this includes positive and negative wavenumbers.
# More precisely, `ks[2:N÷2]` contains the positive wavenumbers, and `ks[N÷2+1:end]`
# contains the corresponding negative wavenumbers.
#
# We now want to compute the wave action spectrum ``n(k) = |\hat{r}(k)|^2 + |\hat{r}(-k)|^2``,
# which is related to the amplitude of the oscillations at the scale ``λ = 2π/k``.

function wave_action_spectrum(ks::AbstractVector, rhat::AbstractVector)
    @assert ks[2] == -ks[end]  # contains positive and negative wavenumbers
    @assert length(ks) == length(rhat)
    N = length(ks)
    if iseven(N)
        Nh = N ÷ 2
        @assert ks[Nh + 1] == -(ks[Nh] + 1)  # wavenumbers change sign after index Nh
    else
        Nh = N ÷ 2 + 1
        @assert ks[Nh + 1] == -ks[Nh]  # wavenumbers change sign after index Nh
    end
    ks_pos = ks[2:Nh]  # only positive wavenumbers
    nk = similar(ks_pos)
    for j ∈ eachindex(ks_pos)
        local k = ks_pos[j]
        i⁺ = 1 + j      # index of coefficient corresponding to wavenumber +k
        i⁻ = N + 1 - j  # index of coefficient corresponding to wavenumber -k
        @assert ks[i⁺] == -ks[i⁻] == k  # verification
        nk[j] = abs2(rhat[i⁺]) + abs2(rhat[i⁻])
    end
    ks_pos, nk
end

ks_pos, nk = wave_action_spectrum(ks, rhat)
nk_normalised = nk ./ ((ϵ * L)^2 / 2)
sum(nk_normalised)  # we expect the sum to be 1

# We can finally plot the final state:

fig = Figure()
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = L"k", ylabel = L"2 \, n(k) / A^2", xlabelsize = 20, ylabelsize = 20)
scatterlines!(ax, ks_pos, nk_normalised)
xlims!(ax, 0.8 * ks_pos[begin], nothing)
ylims!(ax, 1e-30, 1e1)
vlines!(ax, ks_pos[m]; linestyle = :dash, color = :orangered)
fig

# We see that the wave action spectrum is strongly peaked at the wavenumber ``k = 2πm/L``
# (dashed vertical line) corresponding to the perturbation mode ``m`` we chose at the
# beginning (the other peaks are spurious and about 6 orders of magnitude smaller than the
# main peak).
# We also see that the sum ``∑_k n(k)`` (which is basically just the value of the main peak
# in this case) is equal to ``A^2/2``, where ``A = ϵL`` is the amplitude of the initial
# perturbation.
#
# The main conclusion is that, when we perturb a single Kelvin wave mode as we did here,
# that original mode is exactly preserved over time (except for negligible spurious
# effects).

# ## Temporal analysis
#
# We can do something similar to analyse the *temporal* oscillations of the filament.
# For example, we can take the same temporal data we analysed before, corresponding to the
# position of a single filament node:

xt = getindex.(X_probe, 1)  # x positions of a single node over time
yt = getindex.(X_probe, 2)  # y positions of a single node over time
zt = getindex.(X_probe, 3)  # z positions of a single node over time
std(zt) / mean(zt)    # ideally, the z positions shouldn't change over time

# Similarly to before, we now write ``r(t) = x(t) + i y(t)`` and perform an FFT:

inds_t = eachindex(times)[begin:end - 1]  # don't consider the last time to make sure the timestep Δt is constant
rt = @views @. xt[inds_t] + im * yt[inds_t]
Nt = length(rt)           # number of time snapshots
Δt = times[2] - times[1]  # timestep
@assert times[begin:end-1] ≈ range(times[begin], times[end-1]; length = Nt)  # check that times are equispaced
rhat = fft(rt)
@. rhat .= rhat ./ Nt  # normalise FFT
ωs = fftfreq(Nt, 2π / Δt)

ωs_pos, nω = wave_action_spectrum(ωs, rhat)
ωs_normalised = ωs_pos ./ ω_kw  # normalise by expected KW frequency

fig = Figure()
ax = Axis(fig[1, 1]; xscale = log10, yscale = log10, xlabel = L"ω / ω_{\text{kw}}", ylabel = L"n(ω)", xlabelsize = 20, ylabelsize = 20)
scatterlines!(ax, ωs_normalised, nω)
xlims!(ax, 0.8 * ωs_normalised[begin], 1.2 * ωs_normalised[end])
vlines!(ax, 1.0; linestyle = :dash, color = :orangered)
fig

# We see that the temporal spectrum is strongly peaked near the analytical Kelvin wave
# frequency (dashed vertical line).
# Note that the trajectory is not perfectly periodic in time (the ending time doesn't
# exactly match the start time), which can explain non-zero values far from the peak.
