module VortexPasta

include("PaddedArrays/PaddedArrays.jl")  # completely independent of other modules
include("CellLists/CellLists.jl")        # requires PaddedArrays only

include("BasicTypes/BasicTypes.jl")

include("Quadratures/Quadratures.jl")
using .Quadratures
export GaussLegendre

include("Filaments/Filaments.jl")
include("BiotSavart/BiotSavart.jl")
include("Timestepping/Timestepping.jl")

end
