export PointVortexProblem, PointVortexSolver

# TODO:
# - document PointVortexSolver, init
# - define show(::IO, ::PointVortexSolver)
# - checkpoints

using ..BiotSavart:
    Vec2, PointVortices, ParamsPointVortex, SignValue, UP, DOWN,
    PointVortexSigns, PointVortexVelocities, PointVortexScalars,
    PointVortexState, PointVortexCache

include("problem.jl")

# Note: we extend the existing TemporalSchemeCache with methods for PointVortexState.
function init_cache(scheme::TemporalScheme, vortices::PointVortexState)
    Nf = nbuf_filaments(scheme)
    Nv = nbuf_velocities(scheme)
    fc = ntuple(_ -> similar(vortices.points), Val(Nf))
    vc = ntuple(_ -> similar(vortices.points), Val(Nv))
    TemporalSchemeCache(scheme, fc, vc)
end

function Base.resize!(cache::TemporalSchemeCache, vortices::PointVortexState)
    (; fc, vc) = cache
    N = length(vortices)
    for f in fc
        resize!(f, N)
    end
    for v in vc
        resize!(v, N)
    end
    cache
end

function update_values_at_point_vortices!(
        fields::NamedTuple, vortices::PointVortexState, t::Real, iter::AbstractSolver,
    )
    (; to) = iter
    # (0) Apply affect_t! function (optional)
    if iter.affect_t! !== default_callback()
        @timeit to "Affect! (fine-grained)" iter.affect_t!(iter, t)
    end
    # (1) Compute velocities and optionally streamfunction values
    @timeit to "Biot-Savart (full)" BiotSavart.compute_on_nodes!(fields, iter.cache_bs, vortices)
    # add_external_fields!(fields, iter, fs, t, iter.to)
    # apply_forcing!(fields, iter, fs, t, iter.to)
    fields
end

# Case where only the velocity is passed (generally used in internal RK substeps).
# Note: `points` only contains vortex coordinates, it doesn't contain their signs.
# But signs haven't changed from those in iter.vortices (at the start of the timestep).
function update_values_at_point_vortices!(
        vL::AbstractVector, points::AbstractVector, t::Real, iter::AbstractSolver,
    )
    fields = (velocity = vL,)
    (; signs) = iter.vortices
    @assert eachindex(signs) == eachindex(points)
    vortices = PointVortexState(points, signs)
    update_values_at_point_vortices!(fields, vortices, t, iter)
    vL
end

# Note: we keep the `fbase` name for compatibility with vortex filaments (3D)
function advect_point_vortices!(points::AbstractVector, vL::AbstractVector, dt; fbase = nothing)
    _advect_point_vortices!(points, fbase, vL, dt)
end

# Called at the end of a timestep (advance points from their previous locations)
function _advect_point_vortices!(points::AbstractVector, fbase::Nothing, vL::AbstractVector, dt)
    Threads.@threads for i in eachindex(points, vL)
        @inbounds points[i] = points[i] + dt * vL[i]
    end
    points
end

# Called from intermediate RK stages (advance points from a base state)
function _advect_point_vortices!(points::AbstractVector, fbase::AbstractVector, vL::AbstractVector, dt)
    Threads.@threads for i in eachindex(points, vL)
        @inbounds points[i] = fbase[i] + dt * vL[i]
    end
    points
end

"""
    PointVortexSolver{T <: AbstractFloat} <: AbstractSolver{T}

Contains the instantaneous state of a point vortex simulation.

Must be constructed using [`init`](@ref).

# Included fields

Some useful fields included in a `PointVortexSolver` are:

- `prob`: associated [`PointVortexProblem`](@ref) (including point vortex parameters);

- `vortices::PointVortexState`: current state of point vortices in the system (positions and signs);

- `quantities`: a `NamedTuple` containing data on point vortices at the current simulation
  timestep. See **Physical quantities** below for more details.

- `time`: a [`TimeInfo`](@ref) object containing information such as the current time and timestep;

- `to`: a `TimerOutput`, which records the time spent on different functions;

- `cache_bs`: the point vortex cache, which contains data from short- and
  long-range computations.

## Physical quantities

The `quantities` field is a `NamedTuple` containing instantaneous data on point vortices.
For example, `quantities.vs` is the self-induced vortex velocity.

For convenience, if one has a `PointVortexSolver` `iter`, the shortcut `iter.vs` is
equivalent to `iter.quantities.vs` (this also applies to all other quantities).

### List of quantities

- `vs`: self-induced vortex velocity due to Biot–Savart's law.

- `vL`: actual point vortex velocity, including `vs` and potentially other contributions.

- `vf`: external forcing velocity if a forcing scheme has been enabled.

- `vdiss`: external dissipation velocity if a dissipation scheme has been enabled.

- `ψs`: self-induced streamfunction on vortices. This is useful for estimating the kinetic
  energy associated to the induced velocity field.

For convenience, if there is no normal fluid or external forcing, then `vL` is defined as an
alias of `vs` (they're the same object).
"""
struct PointVortexSolver{
        T <: AbstractFloat,
        Problem <: PointVortexProblem{T},
        Vortices <: PointVortexState{T},
        Quantities <: NamedTuple,
        Adaptivity <: AdaptivityCriterion,
        CacheBS <: PointVortexCache,
        CacheTimestepper <: TemporalSchemeCache,
        Timer <: TimerOutput,
        Affect <: Function,
        AffectTime <: Function,
        Callback <: Function,
        AdvectFunction <: Function,
        RHSFunction <: Function,
    } <: AbstractSolver{T}
    prob::Problem
    vortices::Vortices
    quantities::Quantities  # stored fields: (vs, ψs) + optional ones (vn, vf, ...)
    time::TimeInfo
    dtmin::T
    adaptivity::Adaptivity
    cache_bs::CacheBS
    cache_timestepper::CacheTimestepper
    fold_periodic::Bool
    affect!::Affect      # signature: affect!(iter)
    affect_t!::AffectTime  # signature: affect_t!(iter, t) where t is the current time
    callback::Callback   # signature: callback(iter)
    step_diagnostics::Int
    step_adaptivity::Int
    to::Timer
    advect!::AdvectFunction  # function for advecting filaments with a known velocity
    rhs!::RHSFunction        # function for estimating filament velocities (and sometimes streamfunction) from their positions
end

function Base.show(io_in::IO, iter::PointVortexSolver)
    io = IOContext(io_in, :prefix => " │  ")  # this may be used by some `show` functions, for example for the forcing
    print(io, "PointVortexSolver with fields:")
    print(io, "\n ├─ prob: ", iter.prob)
    print(io, "\n ├─ vortices: ", iter.vortices)
    print(io, "\n ├─ vL: ", summary(iter.vL))
    print(io, "\n ├─ vs: ", summary(iter.vs))
    print(io, "\n ├─ ψs: ", summary(iter.ψs))
    print(io, "\n ├─ time: ")
    summary(io, iter.time)
    print(io, "\n ├─ dtmin: ", iter.dtmin)
    print(io, "\n ├─ adaptivity: ", iter.adaptivity)
    print(io, "\n ├─ fold_periodic: ", iter.fold_periodic)
    print(io, "\n ├─ cache_bs: ")
    summary(io, iter.cache_bs)
    print(io, "\n ├─ cache_timestepper: ")
    summary(io, iter.cache_timestepper)
    _maybe_print_function(io, "\n ├─ affect!:", iter.affect!)
    _maybe_print_function(io, "\n ├─ affect_t!:", iter.affect_t!)
    _maybe_print_function(io, "\n ├─ callback:", iter.callback)
    print(io, "\n ├─ step_diagnostics: ", iter.step_diagnostics)
    print(io, "\n ├─ step_adaptivity: ", iter.step_adaptivity)
    # if iter.forcing !== NoForcing()
    #     print(io, "\n ├─ forcing: ", iter.forcing)  # note: this draws a "subtree"
    # end
    # if iter.dissipation !== NoDissipation()
    #     print(io, "\n ├─ dissipation: ", iter.dissipation)  # note: this draws a "subtree"
    # end
    print(io, "\n ├─ advect!: Function")
    print(io, "\n ├─ rhs!: Function")
    print(io, "\n └─ to: ")
    summary(io, iter.to)
end

@doc raw"""
    init(prob::PointVortexProblem, scheme::TemporalScheme; dt::Real, kws...) -> PointVortexProblem

Initialise point vortex problem.

Returns a [`PointVortexSolver`](@ref) which can be advanced in time using
either [`step!`](@ref) or [`solve!`](@ref).

# Mandatory positional arguments

- `prob`: a [`PointVortexProblem`](@ref) containing the problem definition.

- `scheme`: a timestepping scheme (see [`TemporalScheme`](@ref)).

# Mandatory keyword arguments

- `dt::Real`: the simulation timestep.

# Optional keyword arguments

- `adaptivity = NoAdaptivity()`: criterion used for adaptively setting the timestep `dt`.
  See [`AdaptivityCriterion`](@ref) for a list of possible criteria.
  Note that not all criteria are compatible with point vortex simulations
  (e.g. [`AdaptBasedOnSegmentLength`](@ref) only makes sense in 3D).

- `dtmin = 0.0`: minimum `dt` for adaptive timestepping. If `dt < dtmin`, the
  solver is stopped with an error.

- `fold_periodic = true`: if `true` (default), vortices will be recentred onto the main unit cell
  when using periodic boundary conditions. It may be convenient to disable this for
  visualisation purposes. This setting doesn't affect the results, besides possible spatial
  translations of the vortices proportional to the domain period.

- `callback`: a function to be called at the end of each timestep. Its signature must be
  `callback(iter::VortexFilamentSolver)`. **Filaments should not be modified** by this
  function. In that case use `affect!` instead. See notes below for more details.

- `affect!`: similar to `callback`, but allows to modify filaments and other quantities at
  the end of each timestep. Its signature must be `affect!(iter::VortexFilamentSolver)`.
  See notes below for more details.

- `affect_t!`: similar to `affect!`, but allows to modify filaments and other quantities
  _within a timestep_ (e.g. within each Runge–Kutta substep). Its signature must be
  `affect_t!(iter::VortexFilamentSolver, time::Real)`. See notes below for more details.

- `step_diagnostics = 1`: every how many timesteps will diagnostics (from the
  [`Diagnostics`](@ref) module) be computed. This is used to avoid computing certain
  quantities (such as `ψs`) which are only needed for diagnostics (e.g. energy). If this option is used,
  one should call [`can_compute_diagnostics`](@ref) before computing any diagnostics.

- `step_adaptivity = 1`: every how many simulation steps will we adapt the time step `Δt`
  (see `adaptivity` option). When using splitting methods (e.g. [`Strang`](@ref)) combined
  with [`AdaptBasedOnVelocity`](@ref), this can have an extra cost since it will force
  computation of velocities at the start of a time step. For that reason, it can make sense
  to set this equal to `step_diagnostics` (which also force velocity computations).

- `timer = TimerOutput("VortexFilament")`: an optional `TimerOutput` for
  recording the time spent on different functions.
"""
function init(
        prob::PointVortexProblem{T}, scheme::TemporalScheme;
        dt::Real,
        dtmin::Real = zero(T),
        adaptivity::AdaptivityCriterion = NoAdaptivity(),
        fold_periodic::Bool = true,
        callback::Callback = default_callback(),  # by default this is an empty function which just returns `nothing`
        affect!::Affect = default_callback(),
        affect_t!::AffectTime = default_callback(),
        step_diagnostics::Int = 1,
        step_adaptivity::Int = 1,
        timer = TimerOutput("PointVortex"),
    ) where {
        Callback <: Function,
        Affect <: Function,
        AffectTime <: Function,
        T,
    }
    (; tspan, state,) = prob

    vortices = copy(prob.vortices)
    vs = similar(vortices.points)     # vector
    ψs = similar(vortices.points, T)  # scalar

    vL = vs  # actual vortex velocity (alised to vs as long as there are no other external velocities)
    vL_buf = similar(vL)  # this is used as a temporary buffer
    quantities = (; vs, ψs, vL, vL_buf)

    cache_bs = BiotSavart.init_cache(prob.p; timer)
    cache_timestepper = init_cache(scheme, vortices)

    # Wrap functions with the timer, so that timings are estimated each time the function is called.
    advect! = timer(advect_point_vortices!, "Advect vortices")
    rhs! = timer(update_values_at_point_vortices!, "Update fields (RHS)")

    callback_ = _maybe_wrap_function_with_timer(callback, timer, "Callback")
    affect_ = _maybe_wrap_function_with_timer(affect!, timer, "Affect!")
    affect_t_ = _maybe_wrap_function_with_timer(affect_t!, timer, "Affect! (fine-grained)")

    adaptivity_ = possibly_add_max_timestep(adaptivity, dt)
    if adaptivity_ !== NoAdaptivity() && !can_change_dt(scheme)
        throw(ArgumentError(lazy"temporal scheme $scheme doesn't support adaptibility; set `adaptivity = NoAdaptivity()` or choose a different scheme"))
    end

    if state === (;)  # empty state, new simulation
        time = TimeInfo(t = first(tspan), dt = dt, dt_prev = dt)
    else  # with state, restarting from checkpoint
        time = state.time
        if adaptivity === NoAdaptivity()
            time.dt = dt  # use constant timestep given as input
        end
    end

    step_diagnostics > 0 || throw(ArgumentError("step_diagnostics should be positive"))
    step_adaptivity > 0 || throw(ArgumentError("step_adaptivity should be positive"))

    iter = PointVortexSolver(
        prob, vortices, quantities, time, dtmin, adaptivity_, cache_bs,
        cache_timestepper, fold_periodic, affect_, affect_t_, callback_,
        step_diagnostics, step_adaptivity, timer, advect!, rhs!,
    )

    # Verify signature of callback and affect functions
    applicable(callback, iter) || throw(ArgumentError("`callback` function should be callable as `f(iter::VortexFilamentSolver)`"))
    applicable(affect!, iter) || throw(ArgumentError("`affect!` function should be callable as `f(iter::VortexFilamentSolver`)"))
    applicable(affect_t!, iter, time.t) || throw(ArgumentError("`affect_t!` function should be callable as `f(iter::VortexFilamentSolver, t::Real)`"))

    status = finalise_step!(iter; force_computation_of_fields = true)  # make sure all fields (velocity, streamfunction) are computed at initial step
    status == SUCCESS || error("reached status = $status at initialisation")

    iter
end

get_vortices(iter::PointVortexSolver) = iter.vortices  # positions + signs (PointVortexState)
get_vortex_positions(iter::PointVortexSolver) = iter.vortices.points  # positions vector

# Called whenever vortex positions have just been initialised or updated.
function finalise_step!(iter::PointVortexSolver; force_computation_of_fields = false)::SimulationStatus
    (; vortices, time, adaptivity, rhs!, to) = iter
    (; vs, vL, ψs,) = iter.quantities

    @assert eachindex(vortices) == eachindex(vL) == eachindex(vs)
    isempty(vortices) && return NO_VORTICES_LEFT

    iter.affect!(iter)

    adapt_this_timestep = iter.nstep % iter.step_adaptivity == 0
    require_v_and_ψ = can_compute_diagnostics(iter) || force_computation_of_fields
    require_v = require_v_and_ψ || requires_full_velocity(scheme(iter)) || (adapt_this_timestep && requires_full_velocity(iter.adaptivity))

    if require_v_and_ψ
        let fields = (; velocity = vL, streamfunction = ψs)
            rhs!(fields, vortices, time.t, iter)
        end
    elseif require_v
        let fields = (; velocity = vL)
            rhs!(fields, vortices, time.t, iter)
        end
    end

    if adapt_this_timestep
        time.dt_prev = time.dt
        time.dt = estimate_timestep(adaptivity, iter)  # estimate dt for next timestep
    end

    iter.callback(iter)

    SUCCESS
end

function after_advection!(iter::PointVortexSolver)
    (; vortices, to) = iter
    (; points) = vortices
    if iter.fold_periodic
        @timeit to "Periodic fold" begin
            Ls = periods(iter.prob.p)  # box size (periodicity)
            Threads.@threads for i in eachindex(points)
                @inbounds points[i] = Filaments.fold_coordinates_periodic(points[i], Ls)
            end
        end
    end
    iter
end

function step!(iter::PointVortexSolver)::SimulationStatus
    (; vL_buf) = iter.quantities
    _step!(iter, iter.vortices.points, vL_buf)
end

include("diagnostics.jl")
