"""
    MultirateScheme <: TemporalScheme

Abstract type defining a multirate scheme.

The idea is to treat different terms of the evolution equations with different timesteps
(and possibly different integration schemes). Concretely, the fast dynamics is integrated
with a smaller timestep than the slowly-evolving motions.
"""
abstract type MultirateScheme <: TemporalScheme end

include("midpoint.jl")
include("MRI-GARK-ERK33a.jl")
