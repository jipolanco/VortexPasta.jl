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
    @includetest "trefoil.jl"
    @includetest "plots.jl"
    @includetest "synthetic_fields.jl"
    @includetest "forcing.jl"
    @includetest "constants.jl"
    @includetest "vector_of_vector.jl"
    @includetest "padded_arrays.jl"
    @includetest "splines.jl"
    @includetest "filaments.jl"
    @includetest "min_distance.jl"
    @includetest "refinement.jl"
    @includetest "hdf5.jl"
    @includetest "ring.jl"
    @includetest "background_vorticity.jl"
    # -- Tests with timestepping start here -- #
    @includetest "inject_filaments.jl"
    @includetest "ring_lia.jl"
    @includetest "ring_energy.jl"
    @includetest "ring_perturbed.jl"
    @includetest "ring_collision.jl"
    @includetest "ring_stretching.jl"
    @includetest "ring_friction.jl"
    @includetest "checkpoint.jl"
    @includetest "remove_small_filaments.jl"
    @includetest "minimal_energy.jl"
    @includetest "links.jl"
    @includetest "infinite_lines.jl"
    @includetest "imex.jl"
    @includetest "leapfrogging.jl"
    @includetest "kelvin_waves.jl"
    @includetest "forced_lines.jl"
    @includetest "reconnections.jl"
end
