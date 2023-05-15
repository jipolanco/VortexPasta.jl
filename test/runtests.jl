using Test
using VortexPasta

@testset "VortexPasta.jl" begin
    include("filaments.jl")
    include("ring.jl")
    include("ring_collision.jl")
    include("trefoil.jl")
    include("infinite_lines.jl")
    include("plots.jl")
end
