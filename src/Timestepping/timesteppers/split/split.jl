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

include("strang.jl")
