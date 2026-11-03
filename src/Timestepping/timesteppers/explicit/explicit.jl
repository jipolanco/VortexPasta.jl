"""
    ExplicitScheme <: TemporalScheme

Abstract type defining an explicit temporal scheme.
"""
abstract type ExplicitScheme <: TemporalScheme end

include("Euler.jl")
include("midpoint.jl")
include("RK4.jl")
include("SSPRK33.jl")
include("DP5.jl")
include("Vern6.jl")
include("Vern8.jl")
include("Vern9.jl")
