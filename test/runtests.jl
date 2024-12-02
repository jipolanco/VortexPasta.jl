using Test
using InteractiveUtils: versioninfo

# Wraps each test file in a separate module, to avoid definition clashes and to make sure
# that each file can also be run as a standalone script.
macro includetest(path::String)
    modname = Symbol("Mod_" * replace(path, '.' => '_'))
    escname = esc(modname)
    ex = quote
        module $escname
            elapsedtime = time_ns()
            $escname.include($path)
            elapsedtime = (time_ns() - elapsedtime) / 1e9
            printstyled(" * Completed $(rpad($path, 20))"; bold = true, color = :light_green)
            printstyled(" $elapsedtime s\n"; color = :light_green)
        end
    end
    ex.head = :toplevel
    ex
end

println()
@show Base.julia_cmd()
println()
versioninfo()
println()

@info "Running tests with $(Threads.nthreads()) threads"

@testset "VortexPasta.jl" begin
    @includetest "plots.jl"
    @includetest "forcing.jl"
    @includetest "constants.jl"
    @includetest "vector_of_vector.jl"
    @includetest "padded_arrays.jl"
    @includetest "splines.jl"
    @includetest "filaments.jl"
    @includetest "refinement.jl"
    @includetest "hdf5.jl"
    @includetest "ring.jl"
    @includetest "ring_lia.jl"
    @includetest "ring_energy.jl"
    @includetest "ring_perturbed.jl"
    @includetest "ring_collision.jl"
    @includetest "ring_stretching.jl"
    @includetest "background_vorticity.jl"
    @includetest "remove_small_filaments.jl"
    @includetest "trefoil.jl"
    @includetest "links.jl"
    @includetest "infinite_lines.jl"
    @includetest "imex.jl"
    @includetest "leapfrogging.jl"
    @includetest "kelvin_waves.jl"
    @includetest "forced_lines.jl"
    @includetest "min_distance.jl"
    @includetest "reconnections.jl"
end
