# Some environment variables which can influence tests:
#
# - `JULIA_ENABLE_JET_KA_TESTS=true`: enables inference tests based on JET.jl of functions
#   that may call KernelAbstractions.jl kernels at some point. Those kernels are by
#   construction type unstable, so those tests will always fail...
#
# - `JULIA_TESTS_VERBOSE=true`: prints more information during tests, including lots of
#   plots using UnicodePlots.jl and timings from TimerOutputs.jl.

using Test
using InteractiveUtils: versioninfo

# Wraps each test file in a separate module, to avoid definition clashes and to make sure
# that each file can also be run as a standalone script.
macro includetest(testname::String, path::String)
    modname = Symbol("Mod_" * replace(path, '.' => '_'))
    escname = esc(modname)
    ex = quote
        module $escname; end  # create empty module where test will be evaluated (to avoid polluting main module)
        let elapsedtime = time_ns()
            @testset $testname begin
                $escname.include($path)  # evaluate test in empty module
            end
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

@testset "VortexPasta.jl" verbose=true begin
    @includetest("Plots", "plots.jl")
    @includetest("Synthetic fields", "synthetic_fields.jl")
    @includetest("Constants", "constants.jl")
    @includetest("VectorOfVector", "vector_of_vector.jl")
    @includetest("Padded arrays", "padded_arrays.jl")
    @includetest("Cell lists", "cell_lists.jl")
    @includetest("Splines", "splines.jl")
    @includetest("Filaments", "filaments.jl")
    @includetest("Minimum distance", "min_distance.jl")
    @includetest("Refinement", "refinement.jl")
    @includetest("FilamentIO", "hdf5.jl")
    @includetest("Vortex ring", "ring.jl")
    @includetest("Background vorticity", "background_vorticity.jl")
    # -- Tests with timestepping start here -- #
    @includetest("Inject filaments", "inject_filaments.jl")
    @includetest("LIA ring", "ring_lia.jl")
    @includetest("Ring energy", "ring_energy.jl")
    @includetest("Perturbed ring", "ring_perturbed.jl")
    @includetest("Ring collision", "ring_collision.jl")
    @includetest("Ring stretching", "ring_stretching.jl")
    @includetest("Ring friction", "ring_friction.jl")
    @includetest("Forcing", "forcing.jl")
    @includetest("Checkpoints", "checkpoint.jl")
    @includetest("Remove small filaments", "remove_small_filaments.jl")
    @includetest("Minimal energy", "minimal_energy.jl")
    @includetest("Trefoil knot", "trefoil.jl")
    @includetest("Links", "links.jl")
    @includetest("Infinite lines", "infinite_lines.jl")
    @includetest("IMEX", "imex.jl")
    @includetest("Leapfrogging rings", "leapfrogging.jl")
    @includetest("Kelvin waves", "kelvin_waves.jl")
    @includetest("Forced lines", "forced_lines.jl")
    @includetest("Reconnections", "reconnections.jl")
end
