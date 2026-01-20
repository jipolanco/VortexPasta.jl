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
from the default behaviour of [`FilamentIO.write_vtkhdf`](@ref). In particular, this means
that filaments may be broken down into multiple elements. One can pass `periods = nothing`
to disable this behaviour.

Keyword arguments are passed to [`FilamentIO.write_vtkhdf`](@ref).

## Saving velocity and other quantities

One can save other quantities for visualisation by passing an optional `f` function (usually
via a `do` block, as in the example below). See [`FilamentIO.write_vtkhdf`](@ref) and the
examples below for more details.

## Examples

Run simulation and save checkpoint at the end:

```julia
prob = VortexFilamentProblem(...)
iter = init(prob, RK4(); ...)
solve!(iter)
filename = "filaments_\$(iter.nstep).vtkhdf"
save_checkpoint(filename, iter) do io
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
        gtop["version"] = string(VortexPasta.version())
        gtop["julia_version"] = string(VERSION)
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
            g["temporal_scheme"] = string(scheme(iter.cache_timestepper))
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

save_checkpoint(filename, iter::VortexFilamentSolver; kwargs...) = save_checkpoint(identity, filename, iter; kwargs...)

struct LoadedCheckpoint{Filaments, Time <: TimeInfo, Stats <: SimulationStats}
    fs    :: Filaments
    time  :: Time
    stats :: Stats
end

"""
    load_checkpoint(filename, T, method::DiscretisationMethod; read_time = true) -> LoadedCheckpoint

Load simulation state previously written by [`save_checkpoint`](@ref).

The resulting checkpoint can be used to construct a [`VortexFilamentProblem`](@ref).

If `read_time = false`, time information will _not_ be read from the checkpoint file.
This can be used to make sure that the new simulation starts at time `t = 0` (and `nstep = 0`).

## Examples

Restart simulation from checkpoint:

```julia
using VortexPasta.BiotSavart
using VortexPasta.Timestepping
p = ParamsBiotSavart(...)
checkpoint = load_checkpoint("filaments_1234.vtkhdf", Float64, QuinticSplineMethod())
tsim = 2.0  # total simulation time
prob = VortexFilamentProblem(checkpoint, tsim, p)
iter = init(prob, RK4(); ...)
solve!(iter)
```
"""
function load_checkpoint(filename, ::Type{T}, method::DiscretisationMethod; read_time = true) where {T <: AbstractFloat}
    time = TimeInfo()  # this initialises all fields to 0 (t = 0, nstep = 0, ...)
    stats = SimulationStats(T)
    fs = FilamentIO.read_vtkhdf(filename, T, method) do io
        gbase = FilamentIO.root(io)  # "/" group in HDF5 file
        haskey(gbase, "VortexPasta") || error(
            lazy"""HDF5 file "$filename" doesn't contain a /VortexPasta group.
            Are you sure it was written using save_checkpoint?
            Note that files written using FilamentIO.write_vtkhdf cannot be loaded as checkpoints since they don't have enough data.
            """
        )
        gtop = HDF5.open_group(gbase, "VortexPasta")
        let g = HDF5.open_group(gtop, "VortexFilamentSolver")
            if read_time
                let g = HDF5.open_group(g, "Time")
                    _read_all_properties!(g, time)
                end
            end
            let g = HDF5.open_group(g, "SimulationStats")
                _read_all_properties!(g, stats)
            end
        end
    end
    LoadedCheckpoint(fs, time, stats)
end

# Create problem from checkpoint
function VortexFilamentProblem(
        checkpoint::LoadedCheckpoint, tsim::Real, p::ParamsBiotSavart,
    )
    (; fs, time, stats,) = checkpoint
    state = (; time, stats,)  # solver state
    tspan = (time.t, time.t + tsim)
    VortexFilamentProblem(fs, tspan, p, state)
end

# Write all properties of object to HDF5 group
function _write_all_properties!(g, obj)
    foreach(propertynames(obj)) do name
        @inline
        g[string(name)] = getproperty(obj, name)
    end
    nothing
end

function _read_all_properties!(g, obj)
    foreach(propertynames(obj)) do name
        @inline
        T = fieldtype(typeof(obj), name)
        val = if haskey(g, string(name))
            read(g, string(name) => T)::T
        else
            @warn "Property not found in checkpoint file: $name" g
            zero(T)  # if the property doesn't exist (e.g. if the file was written with an old VP version where it didn't exist)
        end
        setproperty!(obj, name, val)
    end
    obj
end
