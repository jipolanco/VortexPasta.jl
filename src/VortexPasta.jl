module VortexPasta

include("BasicTypes/BasicTypes.jl")

include("Quadratures/Quadratures.jl")
using .Quadratures
export GaussLegendreQuadrature

include("Filaments/Filaments.jl")
include("BiotSavart/BiotSavart.jl")

end
