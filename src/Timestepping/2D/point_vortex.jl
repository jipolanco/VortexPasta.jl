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
    PointVortexSolver{T <: AbstractFloat} <: AbstractSolver
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
