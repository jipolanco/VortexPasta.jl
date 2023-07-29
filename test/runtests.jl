using Test
using VortexPasta

@testset "VortexPasta.jl" begin
    include("vector_of_vector.jl")
    include("filaments.jl")
    include("refinement.jl")
    include("hdf5.jl")
    include("ring.jl")
    include("ring_collision.jl")
    include("trefoil.jl")
    include("infinite_lines.jl")
    include("leapfrogging.jl")
    include("kelvin_waves.jl")
    include("min_distance.jl")
    include("reconnections.jl")
    include("plots.jl")
end
