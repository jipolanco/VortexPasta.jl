using Test

# Wraps each test file in a separate module, to avoid definition clashes and to make sure
# that each file can also be run as a standalone script.
macro includetest(path)
    modname = Symbol("Mod_" * replace(path, '.' => '_'))
    ex = quote
        module $(esc(modname))
            $(esc(modname)).include($path)
        end
        using .$(modname)
    end
    ex.head = :toplevel
    ex
end

@testset "VortexPasta.jl" begin
    @includetest "vector_of_vector.jl"
    @includetest "filaments.jl"
    @includetest "refinement.jl"
    @includetest "hdf5.jl"
    @includetest "ring.jl"
    @includetest "ring_collision.jl"
    @includetest "trefoil.jl"
    @includetest "infinite_lines.jl"
    @includetest "leapfrogging.jl"
    @includetest "kelvin_waves.jl"
    @includetest "min_distance.jl"
    @includetest "reconnections.jl"
    @includetest "plots.jl"
end
