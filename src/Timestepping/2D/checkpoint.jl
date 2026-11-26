function save_checkpoint(f::F, filename, iter::PointVortexSolver; kwargs...) where {F <: Function}
    periods = iter.prob.p.Ls  # this currently does nothing in the 2D version of write_vtkhdf
    FilamentIO.write_vtkhdf(filename, iter.vortices; periods, kwargs...) do io
        f(io)  # do optional stuff in /VTKHDF group, e.g. writing other fields such as velocity
        gbase = FilamentIO.root(io)  # "/" group in HDF5 file
        gtop = HDF5.create_group(gbase, "VortexPasta")
        let g = HDF5.create_group(gtop, "ParamsPointVortex")
            BiotSavart.to_hdf5(g, iter.prob.p)
        end
        let g = HDF5.create_group(gtop, "PointVortexSolver")
            let g = HDF5.create_group(g, "Time")
                g["simulation_time_span"] = collect(iter.prob.tspan)
                g["simulation_dtmin"] = iter.dtmin
                _write_all_properties!(g, iter.time)
            end
            # let g = HDF5.create_group(g, "SimulationStats")
            #     _write_all_properties!(g, iter.stats)
            # end
            g["temporal_scheme"] = string(scheme(iter.cache_timestepper))
            g["adaptivity"] = string(iter.adaptivity)
            g["annihilate"] = string(Reconnections.criterion(iter.annihilate))
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
            g["forcing_type"] = string(typeof(iter.forcing))
            g["dissipation_type"] = string(typeof(iter.dissipation))
        end
    end
end

save_checkpoint(filename, iter::PointVortexSolver; kwargs...) = save_checkpoint(identity, filename, iter; kwargs...)

struct LoadedCheckpoint2D{Vortices, Time <: TimeInfo}
    vortices::Vortices
    time::Time
end

function load_checkpoint(func::F, filename, ::Type{T}; read_time = true) where {F <: Function, T <: AbstractFloat}
    time = TimeInfo()  # this initialises all fields to 0 (t = 0, nstep = 0, ...)
    # stats = VortexFilamentStats(T)
    vortices = FilamentIO.read_vtkhdf(filename, T) do io
        gbase = FilamentIO.root(io)  # "/" group in HDF5 file
        haskey(gbase, "VortexPasta") || error(
            lazy"""HDF5 file '$filename' doesn't contain a /VortexPasta group.
            Are you sure it was written using save_checkpoint?
            Note that files written using FilamentIO.write_vtkhdf cannot be loaded as checkpoints since they don't have enough data.
            """
        )
        gtop = HDF5.open_group(gbase, "VortexPasta")
        haskey(gtop, "PointVortexSolver") || error(
            lazy"""HDF5 file '$filename' doesn't contain a /VortexPasta/PointVortexSolver group.
            Are you sure it comes from a 2D point vortex simulation?
            """
        )
        let g = HDF5.open_group(gtop, "PointVortexSolver")
            if read_time
                let g = HDF5.open_group(g, "Time")
                    _read_all_properties!(g, time)
                end
            end
            # let g = HDF5.open_group(g, "SimulationStats")
            #     _read_all_properties!(g, stats)
            # end
        end
        func(io)
    end
    LoadedCheckpoint2D(vortices, time)
end

load_checkpoint(filename, ::Type{T}; kws...) where {T <: AbstractFloat} =
    load_checkpoint(Returns(nothing), filename, T; kws...)

# Create problem from checkpoint
function PointVortexProblem(
        checkpoint::LoadedCheckpoint2D, tsim::Real, p::ParamsPointVortex,
    )
    (; vortices, time) = checkpoint
    state = (; time)  # solver state
    tspan = (time.t, time.t + tsim)
    PointVortexProblem(vortices, tspan, p, state)
end
