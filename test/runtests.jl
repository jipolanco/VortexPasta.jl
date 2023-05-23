using Test
using VortexPasta

@testset "VortexPasta.jl" begin
    include("filaments.jl")
    include("ring.jl")
    include("ring_collision.jl")
    include("trefoil.jl")
    include("infinite_lines.jl")
    include("kelvin_waves.jl")
    include("plots.jl")
end
