export PointVortexProblem, PointVortexSolver

# TODO:
# - document PointVortexSolver, init
# - define show(::IO, ::PointVortexSolver)

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

# Note: we keep the `fbase` name for compatibility with vortex filaments (3D)
function advect_point_vortices!(points::AbstractVector, vL::AbstractVector, dt; fbase = points)
    Threads.@threads for i in eachindex(points, vL, fbase)
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
    } <: AbstractSolver
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
    (; Ls,) = prob.p

    vortices = copy(prob.vortices)
    vs = similar(vortices.points)     # vector
    ψs = similar(vortices.points, T)  # scalar
    quantities = (; vs, ψs,)

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

    # TODO: implement
    # status = finalise_step!(iter; force_computation_of_fields = true)  # make sure all fields (velocity, streamfunction) are computed at initial step
    # status == SUCCESS || error("reached status = $status at initialisation")

    iter
end
