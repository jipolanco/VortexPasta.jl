using Test
using VortexFilamentEwald

@testset "VortexFilamentEwald.jl" begin
    include("filaments.jl")
    include("trefoil.jl")
    include("plots.jl")
end
