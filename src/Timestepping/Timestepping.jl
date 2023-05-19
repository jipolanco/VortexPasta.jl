"""
    Timestepping

Module defining timestepping solvers for vortex filament simulations.
"""
module Timestepping

export init, solve!, step!, VortexFilamentProblem, RK4

using ..BasicTypes: Vec3
using ..Filaments: Filaments, AbstractFilament, nodes, RefinementCriterion, BasedOnCurvature
using ..BiotSavart:
    BiotSavart,
    ParamsBiotSavart,
    BiotSavartCache,
    VectorOfFilaments,
    VectorOfPositions,
    VectorOfVelocities,
    AllFilamentVelocities

# Reuse same init and solve! functions from the SciML ecosystem, to avoid clashes.
# See https://docs.sciml.ai/CommonSolve/stable/
import CommonSolve: init, solve!
using RecursiveArrayTools: VectorOfArray  # for convenience, to deal with filament velocities

using TimerOutputs: TimerOutput, @timeit

abstract type AbstractProblem end
abstract type AbstractSolver end

include("timesteppers/timesteppers.jl")

struct VortexFilamentProblem{Filaments <: VectorOfFilaments} <: AbstractProblem
    fs    :: Filaments
    tspan :: NTuple{2, Float64}
    p     :: ParamsBiotSavart
end

mutable struct VortexFilamentSolver{
        Problem <: VortexFilamentProblem,
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfArray{<:Vec3},
        Refinement <: RefinementCriterion,
        CacheBS <: BiotSavartCache,
        CacheTimestepper <: TemporalSchemeCache,
        Timer <: TimerOutput,
    } <: AbstractSolver
    const prob :: Problem
    const fs   :: Filaments
    const vs   :: Velocities
    t    :: Float64
    dt   :: Float64
    const refinement          :: Refinement
    const cache_bs            :: CacheBS
    const cache_timestepper   :: CacheTimestepper
    const to :: Timer
end

function init(
        prob::VortexFilamentProblem, scheme::ExplicitTemporalScheme;
        alias_u0 = true,  # same as in OrdinaryDiffEq.jl
        dt,
        refinement = BasedOnCurvature(0.35),
        timer = TimerOutput(),
    )
    (; fs, tspan,) = prob
    vs_data = [similar(nodes(f)) for f ∈ fs] :: AllFilamentVelocities
    vs = VectorOfArray(vs_data)
    fs_sol = alias_u0 ? fs : copy.(fs)
    cache_bs = BiotSavart.init_cache(prob.p; timer)
    cache_timestepper = init_cache(scheme, fs, vs)
    t = first(tspan)
    iter = VortexFilamentSolver(
        prob, fs_sol, vs, t, dt, refinement,
        cache_bs, cache_timestepper, timer,
    )
    iter
end

function vortex_velocities!(
        vs::VectorOfArray,
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    (; cache_bs,) = iter
    BiotSavart.velocity_on_nodes!(vs.u, cache_bs, fs)
    vs
end

function vortex_velocities!(
        fs::VectorOfFilaments,
        iter::VortexFilamentSolver,
    )
    vortex_velocities!(iter.vs, fs, iter)
end

function _advect_filament!(
        f::AbstractFilament, vs::VectorOfVelocities, dt::Real;
        fbase = f, L_fold = nothing, refinement = nothing,
    )
    Xs = nodes(f)
    Xs_base = nodes(fbase)
    @assert eachindex(Xs) == eachindex(Xs_base) == eachindex(vs)
    @inbounds for i ∈ eachindex(Xs)
        Xs[i] = Xs_base[i] + dt * vs[i]
    end
    if L_fold !== nothing
        Filaments.fold_periodic!(f, L_fold)
    end
    if refinement !== nothing
        nref = Filaments.refine!(f, refinement)
        if nref !== (0, 0)
            # @info "Added/removed $nref nodes"
            if sum(nref) ≠ 0
                resize!(vs, length(f))
            end
        end
    end
    Filaments.update_coefficients!(f)
    f
end

function _advect_filaments!(fs, vs, dt; fbase = fs, kws...)
    for (f, f₀, v) ∈ zip(fs, fbase, vs)
        _advect_filament!(f, v, dt; fbase = f₀, kws...)
    end
    fs
end

function solve!(iter::VortexFilamentSolver)
    t_end = iter.prob.tspan[2]
    while iter.t < t_end
        step!(iter)
    end
    iter
end

function step!(iter::VortexFilamentSolver)
    (; fs, dt, prob, refinement,) = iter
    vs = _update_velocities!(vortex_velocities!, _advect_filaments!, iter.cache_timestepper, iter)
    L_fold = prob.p.common.Ls  # box size (periodicity)
    _advect_filaments!(fs, vs, dt; L_fold, refinement)
    iter.t += dt
    iter
end

end
