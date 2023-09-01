"""
    ExplicitScheme

Abstract type defining an explicit temporal scheme.
"""
abstract type ExplicitScheme <: TemporalScheme end

include("Euler.jl")
include("RK4.jl")
include("SSPRK33.jl")
include("DP5.jl")
