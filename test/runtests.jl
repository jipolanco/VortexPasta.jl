using Test
using VortexPasta

@testset "VortexPasta.jl" begin
    include("filaments.jl")
    include("ring.jl")
    include("trefoil.jl")
    include("plots.jl")
end
