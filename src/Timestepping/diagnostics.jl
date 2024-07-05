# Overload diagnostics functions for convenience, so that we can simply pass the current
# state of the solver to get diagnostics.

using ..Diagnostics: Diagnostics
using ..Filaments: Filaments

"""
    Diagnostics.kinetic_energy_from_streamfunction(iter::VortexFilamentSolver; quad = nothing)

Estimate kinetic energy from current simulation state.

See [`Diagnostics.kinetic_energy_from_streamfunction`](@ref) for details.
"""
function Diagnostics.kinetic_energy_from_streamfunction(iter::VortexFilamentSolver; kws...)
    (; ψs, fs, external_forcing, t,) = iter
    Ls = BiotSavart.periods(iter.prob.p)
    Γ = BiotSavart.circulation(iter.prob.p)
    E = Diagnostics.kinetic_energy_from_streamfunction(fs, ψs, Γ, Ls; kws...)
    # Add kinetic energy of external velocity field, if available.
    # Note that we only do this if we also included the streamfunction, since otherwise
    # we don't have enough information to estimate the total kinetic energy.
    if external_forcing.velocity !== nothing && external_forcing.streamfunction !== nothing
        E += Diagnostics.kinetic_energy_of_periodic_velocity_field(Ls) do x⃗
            external_forcing.velocity(x⃗, t)
        end
    end
    E
end

"""
    Diagnostics.kinetic_energy_nonperiodic(iter::VortexFilamentSolver; quad = nothing)

Estimate kinetic energy from current simulation state.

See [`Diagnostics.kinetic_energy_nonperiodic`](@ref) for details.

!!! warning

    This function uses an energy definition that is invalid in periodic domains.
    Moreover, even in non-periodic domains, it may display artificial energy fluctuations in
    cases where energy should be conserved. Therefore, it is recommended to **always**
    use [`kinetic_energy_from_streamfunction`](@ref Diagnostics.kinetic_energy_from_streamfunction(::VortexFilamentSolver)), which doesn't have those issues.

"""
function Diagnostics.kinetic_energy_nonperiodic(iter::VortexFilamentSolver; kws...)
    (; vs, fs,) = iter
    Ls = BiotSavart.periods(iter.prob.p)
    BiotSavart.domain_is_periodic(iter.prob.p) &&
        @warn(lazy"`kinetic_energy_nonperiodic` should only be called when working with non-periodic domains (got Ls = $Ls)")
    Γ = BiotSavart.circulation(iter.prob.p)
    Diagnostics.kinetic_energy_nonperiodic(fs, vs, Γ; kws...)
end

# Note: filament_length is actually defined in the Filaments module, but we document its
# alias in Diagnostics just for consistency with the other diagnostics (it doesn't make any
# difference really!).
"""
    Diagnostics.filament_length(iter::VortexFilamentSolver; quad = nothing) -> Real

Estimate total length of all filaments in a simulation.

See [`Diagnostics.filament_length`](@ref) for details.
"""
function Diagnostics.filament_length(iter::VortexFilamentSolver; kws...)
    Diagnostics.filament_length(iter.fs; kws...)
end

"""
    Diagnostics.vortex_impulse(iter::VortexFilamentSolver; quad = nothing) -> Vec3

Estimate total normalised impulse of all filaments in a simulation.

See [`Diagnostics.vortex_impulse`](@ref) for details.
"""
function Diagnostics.vortex_impulse(iter::VortexFilamentSolver; kws...)
    Diagnostics.vortex_impulse(iter.fs; kws...)
end

"""
    Diagnostics.helicity(iter::VortexFilamentSolver; quad = nothing) -> Real

Compute helicity of the instantaneous vortex configuration in a simulation.

See [`Diagnostics.helicity`](@ref) for details.
"""
function Diagnostics.helicity(iter::VortexFilamentSolver; kws...)
    Γ = BiotSavart.circulation(iter.prob.p)
    Diagnostics.helicity(iter.fs, iter.vs, Γ; kws...)
end

"""
    Diagnostics.stretching_rate(iter::VortexFilamentSolver; quad = nothing) -> Real

Compute stretching rate of the instantaneous vortex configuration in a simulation.

This corresponds to the instantaneous rate of increase (or decrease) of total vortex length
in the simulation. It has units of a velocity (``L T^{-1}``).

See [`Diagnostics.stretching_rate`](@ref) for details.
"""
function Diagnostics.stretching_rate(iter::VortexFilamentSolver; kws...)
    (; fs, vs,) = iter
    Diagnostics.stretching_rate(fs, vs; kws...)
end

# This allows passing a VortexFilamentSolver to energy_spectrum / energy_spectrum!.
Diagnostics.get_long_range_cache(iter::VortexFilamentSolver) = iter.cache_bs.longrange
