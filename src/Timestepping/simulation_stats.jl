"""
    SimulationStats{T <: AbstractFloat}

Contains accumulated statistics of different events occurring since the beginning of a
simulation.

Note: for performance reasons, all of the lengths below are estimated using a straight
segment approximation.

# Available fields

- `reconnection_count::Int`: **total number of reconnections**;

- `reconnection_length_loss::T`: accumulated **decrease of filament length due to reconnections**.
  This is estimated from each reconnection event as the difference between the local vortex
  lengths before and after the reconnection, using the straight segment approximation.

- `reconnection_passes::Int`: **total number of reconnection passes**. This is the total
  number of reconnection scans since the beginning of the simulation. Note that the
  [`ReconnectBasedOnDistance`](@ref) and [`ReconnectFast`](@ref) accept a `max_passes` option,
  thus allowing multiple reconnection passes per timestep;

- `filaments_removed_count::Int`: **total number of removed filaments** (this generally happens
  when the number of discretisation points becomes too small for the spatial discretisation
  to work);

- `filaments_removed_length::T`: **total length of removed filaments**. Note that this length is
  estimated using a straight segment approximation (no quadratures). This is because
  filaments are removed when they can no longer be represented using a continuous
  interpolation function;

- `refinement_length_loss::T`: **loss of filament length due to (de)refinement**. When
  performing refinement, filament nodes may be removed to respect the minimal allowed
  distance between nodes. This variable contains the resulting accumulated loss of vortex
  length. In principle this quantity can also be negative, which could happen if nodes are only
  inserted and not removed.
"""
mutable struct SimulationStats{T <: AbstractFloat}
    reconnection_count       :: Int
    reconnection_length_loss :: T
    reconnection_passes      :: Int
    filaments_removed_count  :: Int
    filaments_removed_length :: T
    refinement_length_loss   :: T
end

SimulationStats(::Type{T}) where {T} = SimulationStats(0, zero(T), 0, 0, zero(T), zero(T))

function Base.show(io::IO, stats::SimulationStats)
    print(io, "SimulationStats:")
    print(io, "\n - reconnection_count        = ", stats.reconnection_count)
    print(io, "\n - reconnection_length_loss  = ", stats.reconnection_length_loss)
    print(io, "\n - reconnection_passes       = ", stats.reconnection_passes)
    print(io, "\n - filaments_removed_count   = ", stats.filaments_removed_count)
    print(io, "\n - filaments_removed_length  = ", stats.filaments_removed_length)
    print(io, "\n - refinement_length_loss    = ", stats.refinement_length_loss)
end

Base.summary(io::IO, stats::SimulationStats) = print(io, typeof(stats), "(...)")
