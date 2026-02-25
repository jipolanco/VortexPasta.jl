module VortexPasta

version() = pkgversion(@__MODULE__)::VersionNumber

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

include("VectorsOfVectors/VectorsOfVectors.jl")
include("SyntheticFields/SyntheticFields.jl")
include("Forcing/Forcing.jl")
include("Timestepping/Timestepping.jl")

# This is used in tests
const ALL_MODULES = (
    PaddedArrays, PredefinedCurves, CellLists, Constants,
    Quadratures, Filaments, FindNearbySegments, FilamentIO,
    Reconnections, BiotSavart, Diagnostics, VectorsOfVectors,
    SyntheticFields, Forcing,
    Timestepping,
)

## Precompilation

# For some reason precompilation fails on Julia 1.11, so we only do it for Julia ≥ v1.12.
using PrecompileTools: @setup_workload, @compile_workload

if VERSION ≥ v"1.12"
    @setup_workload begin
        using .PredefinedCurves
        using .Filaments
        using .BiotSavart
        using .Reconnections
        using .Timestepping
        using .Diagnostics

        get_quadrature(::CubicSplineMethod) = GaussLegendre(2)
        get_quadrature(::QuinticSplineMethod) = GaussLegendre(3)

        function run_simulation(; splitting, method)
            quadrature = get_quadrature(method)
            params = ParamsBiotSavart(;
                Γ = 1.2, a = 1.0e-6, Δ = 1 / 4,
                splitting,
                backend_short = CellListsBackend(2),
                backend_long = NonuniformFFTsBackend(),
                quadrature,
                quadrature_near_singularity = quadrature,
                lia_segment_fraction = 0.1,
            )

            # Initialise vortex ring
            S = define_curve(Ring(); scale = π / 3)
            f = Filaments.init(S, ClosedFilament, 32, method)
            fs = [f]
            l_min = minimum_node_distance(fs)

            # Initialise and run simulation
            tspan = (0.0, 0.01)
            prob = VortexFilamentProblem(fs, tspan, params)
            energy = Float64[]

            function callback(iter)
                E = Diagnostics.kinetic_energy(iter)
                push!(energy, E)
                nothing
            end

            iter = init(
                prob, RK4();
                dt = 0.001,
                adaptivity = AdaptBasedOnSegmentLength(0.5),
                refinement = RefineBasedOnSegmentLength(0.75 * l_min),
                reconnect = ReconnectFast(l_min),
                callback,
            )

            step!(iter)  # perform a single timestep
            ks, Ek = Diagnostics.energy_spectrum(iter)

            nothing
        end

        Ls = (1, 1, 1) .* 2π
        Ns = (1, 1, 1) .* 32
        splittings = (
            GaussianSplitting(; Ls, Ns, β = 2.0),
            KaiserBesselSplitting(; Ls, Ns, β = 8.0),
        )
        methods = (
            CubicSplineMethod(), QuinticSplineMethod(),
        )

        for splitting in splittings, method in methods
            @compile_workload run_simulation(; splitting, method)
        end
    end
end

end
