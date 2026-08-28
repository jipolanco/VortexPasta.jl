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

include("strang.jl")
include("strang4.jl")

include("hasimoto.jl")
include("strang4_hasimoto.jl")
