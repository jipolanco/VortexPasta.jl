using Test
using VortexFilamentEwald

@testset "VortexFilamentEwald.jl" begin
    include("filaments.jl")
    include("plots.jl")
end
