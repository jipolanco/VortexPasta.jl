module VortexPasta

include("PaddedArrays/PaddedArrays.jl")  # completely independent of other modules
include("PredefinedCurves/PredefinedCurves.jl")  # completely independent of other modules

include("CellLists/CellLists.jl")        # requires PaddedArrays only

include("Quadratures/Quadratures.jl")
using .Quadratures
export GaussLegendre, AdaptiveTanhSinh

include("Constants/Constants.jl")

include("Filaments/Filaments.jl")
include("FindNearbySegments/FindNearbySegments.jl")
include("FilamentIO/FilamentIO.jl")

include("Reconnections/Reconnections.jl")

include("BiotSavart/BiotSavart.jl")
include("Diagnostics/Diagnostics.jl")

include("Containers/Containers.jl")  # VectorOfVectors
include("Forcing/Forcing.jl")
include("Timestepping/Timestepping.jl")

# This is used in tests
const ALL_MODULES = (
    PaddedArrays, PredefinedCurves, CellLists, Constants,
    Quadratures, Filaments, FindNearbySegments, FilamentIO,
    Reconnections, BiotSavart, Diagnostics, Containers, Timestepping,
)

## Precompilation
using PrecompileTools: @setup_workload, @compile_workload

@setup_workload begin
    using .PredefinedCurves
    using .Filaments
    using .BiotSavart
    using .Reconnections
    using .Timestepping
    using .Diagnostics

    ## Set-up a vortex filament simulation
    # Grid-related parameters
    Ls = (1, 1, 1) .* 2π
    Ns = (1, 1, 1) .* 32
    kmax = minimum(splat((N, L) -> (N ÷ 2) * 2π / L), zip(Ns, Ls))
    α = kmax / 5
    rcut = 4 * sqrt(2) / α

    # Physical vortex parameters
    Γ = 1.2
    a = 1e-6
    Δ = 1/4  # full core

    function callback(iter)
        E = Diagnostics.kinetic_energy_from_streamfunction(iter)
        nothing
    end

    @compile_workload begin
        params_bs = ParamsBiotSavart(;
            Γ, a, Δ,
            α, rcut, Ls, Ns,
            backend_short = CellListsBackend(2),
            backend_long = NonuniformFFTsBackend(),
            quadrature = GaussLegendre(2),
        )

        # Initialise vortex ring
        S = define_curve(Ring(); scale = π / 3)
        f = Filaments.init(S, ClosedFilament, 32, CubicSplineMethod())
        fs = [f]
        l_min = minimum_node_distance(fs)

        # Initialise and run simulation
        tspan = (0.0, 0.01)
        prob = VortexFilamentProblem(fs, tspan, params_bs)
        iter = init(
            prob, RK4();
            dt = 0.001,
            adaptivity = AdaptBasedOnSegmentLength(1.0),
            refinement = RefineBasedOnSegmentLength(0.75 * l_min),
            reconnect = ReconnectBasedOnDistance(l_min),
            callback,
        )
        solve!(iter)
    end
end

end
