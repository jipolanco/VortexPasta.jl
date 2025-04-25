using HDF5: HDF5

"""
    save_checkpoint([f::Function], filename, iter::VortexFilamentSolver; kwargs...)

Save checkpoint to VTKHDF file (usually with `.vtkhdf` extension).

A checkpoint file is a VTKHDF file containing the current state of the solver.
This includes things like vortex locations, the current time and timestep, and solver parameters.
The file can be opened using ParaView allowing to visualise the current state (vortex
locations and optionally other quantities).

In fact it is a VTKHDF file saved using [`FilamentIO.write_vtkhdf`](@ref) which contains
extra information allowing to restart a simulation from a checkpoint while preserving
continuity as much as possible.

To ease visualisations, filaments are folded onto the periodic box by default, which differs
from the default behaviour of [`Filaments.write_vtkhdf`](@ref). In particular, this means
that filaments may be broken down into multiple elements. One can pass `periods = nothing`
to disable this behaviour.

Keyword arguments are passed to [`FilamentIO.write_vtkhdf`](@ref).

## Saving velocity and other quantities

One can save other quantities for visualisation by passing an optional `f` function (usually
via a `do` block, as in the example below). See [`FilamentIO.write_vtkhdf`](@ref) for more details.

## Examples

Run simulation and save checkpoint at the end:

```julia
prob = VortexFilamentProblem(...)
iter = init(prob, RK4(); ...)
solve!(iter)
save_checkpoint("checkpoint.vtkhdf", iter) do io
    io["velocity_s"] = iter.vs  # save superfluid velocity on vortices for visualisation
end
```
"""
function save_checkpoint(f::F, filename, iter::VortexFilamentSolver; kwargs...) where {F <: Function}
    periods = iter.prob.p.Ls  # for convenience, by default we fold vortex locations onto periodic box
    FilamentIO.write_vtkhdf(filename, iter.fs; periods, kwargs...) do io
        f(io)  # do optional stuff in /VTKHDF group, e.g. writing other fields such as velocity
        gbase = FilamentIO.root(io)  # "/" group in HDF5 file
        gtop = HDF5.create_group(gbase, "VortexPasta")
        let g = HDF5.create_group(gtop, "BiotSavartParams")
            BiotSavart.to_hdf5(g, iter.prob.p)
        end
        let g = HDF5.create_group(gtop, "VortexFilamentSolver")
            let g = HDF5.create_group(g, "Time")
                g["simulation_time_span"] = collect(iter.prob.tspan)
                g["simulation_dtmin"] = iter.dtmin
                _write_all_properties!(g, iter.time)
            end
            let g = HDF5.create_group(g, "SimulationStats")
                _write_all_properties!(g, iter.stats)
            end
            g["refinement"] = string(iter.refinement)
            g["adaptivity"] = string(iter.adaptivity)
            g["reconnect"] = string(Reconnections.criterion(iter.reconnect))
            g["fast_term"] = string(iter.fast_term)
            g["LIA"] = iter.LIA
            g["fold_periodic"] = iter.fold_periodic
            # Write TimerOutputs timer replacing all unicode characters with ASCII ones
            # since h5dump/h5ls don't correctly print unicode (as of HDF5 1.14.6).
            g["timer"] = let io = IOBuffer()
                show(io, iter.to; linechars = :ascii)
                s = String(take!(io))
                replace(s, 'μ' => 'u', '–' => '-')
            end
            g["with_affect"] = iter.affect! !== identity
            g["with_callback"] = iter.callback !== identity
            g["with_external_velocity"] = iter.external_fields.velocity !== nothing
            g["with_external_streamfunction"] = iter.external_fields.streamfunction !== nothing
            g["with_stretching_velocity"] = iter.stretching_velocity !== nothing
            g["forcing_type"] = string(typeof(iter.forcing))
        end
    end
    nothing
end

# Write all properties of object to HDF5 group
function _write_all_properties!(g, obj)
    foreach(propertynames(obj)) do name
        @inline
        g[string(name)] = getproperty(obj, name)
    end
    nothing
end

save_checkpoint(filename, iter::VortexFilamentSolver; kwargs...) = save_checkpoint(identity, filename, iter; kwargs...)
