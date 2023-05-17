"""
    Timestepping

Module defining timestepping solvers for vortex filament simulations.
"""
module Timestepping

export init, solve!

using ..Filaments: Filaments, AbstractFilament, nodes, RefinementCriterion, BasedOnCurvature

using ..BiotSavart:
    BiotSavart,
    ParamsBiotSavart,
    BiotSavartCache,
    init_cache,
    velocity_on_nodes!,
    VectorOfFilaments,
    VectorOfPositions,
    VectorOfVelocities

include("timesteppers/timesteppers.jl")

# Reuse same init and solve! functions from the SciML ecosystem, to avoid clashes.
# See https://docs.sciml.ai/CommonSolve/stable/
import CommonSolve: init, solve!

# function vortex_velocity!(vs::VectorOfVelocities, fs::VectorOfFilaments)
# end

struct VortexFilamentProblem{Filaments <: VectorOfFilaments}
    fs    :: Filaments
    tspan :: NTuple{2, Float64}
    ps    :: ParamsBiotSavart
end

mutable struct VortexFilamentSolver{
        Problem <: VortexFilamentProblem,
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfVelocities,
        Refinement <: RefinementCriterion,
        CacheBS <: BiotSavartCache,
    }
    const prob :: Problem
    fs   :: Filaments
    vs   :: Velocities
    t    :: Float64
    dt   :: Float64
    const refinement :: Refinement
    const cache_bs :: CacheBS
end

function init(
        prob::VortexFilamentProblem, scheme::ExplicitTemporalScheme;
        alias_u0 = true,  # same as in OrdinaryDiffEq.jl
        dt,
        refinement = BasedOnCurvature(0.35),
    )
    (; fs, tspan,) = prob
    vs = [similar(nodes(f)) for f ∈ fs] :: VectorOfVelocities
    fs_sol = alias_u0 ? fs : copy.(fs)
    biot_savart_cache = init_cache(prob.ps)
    iter = VortexFilamentSolver(
        prob, fs_sol, vs, tspan[1], dt, refinement,
        biot_savart_cache,
    )
    iter
end

function solve!(solver::VortexFilamentSolver)
    solver
end

end
