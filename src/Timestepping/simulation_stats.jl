"""
    SimulationStats{T <: AbstractFloat}

Contains accumulated statistics of different events occurring since the beginning of a
simulation.

# Available fields

- `reconnection_count::Int`: **total number of reconnections**;

- `reconnection_length_loss::T`: accumulated **decrease of filament length due to reconnections**.
  This is estimated from each reconnection event as the difference between the local vortex
  lengths before and after the reconnection, using the straight segment approximation.

- `filaments_removed_count::Int`: **total number of removed filaments** (this generally happens
  when the number of discretisation points becomes too small for the spatial discretisation
  to work);

- `filaments_removed_length::T`: **total length of removed filaments**. Note that this length is
  estimated using a straight segment approximation (no quadratures). This is because
  filaments are removed when they can no longer be represented using a continuous
  interpolation function.
"""
mutable struct SimulationStats{T <: AbstractFloat}
    reconnection_count       :: Int
    reconnection_length_loss :: T
    filaments_removed_count  :: Int
    filaments_removed_length :: T
end

SimulationStats(::Type{T}) where {T} = SimulationStats(0, zero(T), 0, zero(T))

function Base.show(io::IO, stats::SimulationStats)
    print(io, "SimulationStats:")
    print(io, "\n - reconnection_count        = ", stats.reconnection_count)
    print(io, "\n - reconnection_length_loss  = ", stats.reconnection_length_loss)
    print(io, "\n - filaments_removed_count   = ", stats.filaments_removed_count)
    print(io, "\n - filaments_removed_length  = ", stats.filaments_removed_length)
end

Base.summary(io::IO, stats::SimulationStats) = print(io, typeof(stats), "(...)")
