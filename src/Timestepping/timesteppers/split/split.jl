"""
    SplittingScheme <: TemporalScheme

Abstract type defining a splitting scheme (such as Strang splitting).

Like IMEX ([`ImplicitExplicitScheme`](@ref)) and multirate ([`MultirateScheme`](@ref))
schemes, the idea is to split the Biot–Savart integral into two terms, respectively
governing the "fast" dynamics (usually the local term) and the "slow" dynamics (non-local
interactions).

The evolution due to both terms is approximated using some kind of Runge–Kutta scheme.
"""
abstract type SplittingScheme <: TemporalScheme end

_check_nsubsteps(fast::TemporalScheme, nsubsteps) = nothing

# Splitting schemes don't need the full velocity at time t.
# That is, the velocity passed to _update_velocities! is ignored, so we don't need to
# precompute the velocity.
requires_full_velocity(::SplittingScheme) = false

function splitting_advance_slow!(slow::TemporalScheme, ftmp, vtmp, iter, τ, rhs_full!::F, advect!::G, cache, cdt) where {F, G}
    rhs!(args...) = rhs_full!(args...; component = Val(:slow))  # function to compute RHS associated to slow term
    dτ = cdt
    rhs!(vtmp, ftmp, τ, iter)  # stage 1 (of e.g. RK4)
    update_velocities!(
        vtmp, rhs!, advect!, cache, iter;
        resize_cache = false, t = τ, dt = dτ, fs = ftmp,
    )
    advect!(ftmp, vtmp, dτ; fbase = ftmp)
    τ + dτ
end

function splitting_advance_fast!(fast::TemporalScheme, ftmp, vtmp, iter, τ, rhs_full!::F, advect!::G, cache, cdt, nsubsteps) where {F, G}
    rhs!(args...) = rhs_full!(args...; component = Val(:fast))  # function to compute RHS associated to fast term
    dτ = cdt / nsubsteps  # timestep in each substep
    for _ ∈ 1:nsubsteps
        rhs!(vtmp, ftmp, τ, iter)
        update_velocities!(
            vtmp, rhs!, advect!, cache, iter;
            resize_cache = false, t = τ, dt = dτ, fs = ftmp,
        )
        advect!(ftmp, vtmp, dτ; fbase = ftmp)
        τ += dτ
    end
    τ
end

include("hasimoto.jl")
include("strang.jl")
include("strang4.jl")
