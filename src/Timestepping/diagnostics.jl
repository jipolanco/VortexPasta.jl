# Overload diagnostics functions for convenience, so that we can simply pass the current
# state of the solver to get diagnostics.

using ..Diagnostics: Diagnostics
using ..Filaments: Filaments

function Diagnostics.kinetic_energy_from_streamfunction(iter::VortexFilamentSolver; kws...)
    (; ψs, fs, external_fields, t,) = iter
    Ls = BiotSavart.periods(iter.prob.p)
    E = Diagnostics.kinetic_energy_from_streamfunction(fs, ψs, iter.prob.p; kws...)
    # Add kinetic energy of external velocity field, if available.
    # Note that we only do this if we also included the streamfunction, since otherwise
    # we don't have enough information to estimate the total kinetic energy.
    if external_fields.velocity !== nothing && external_fields.streamfunction !== nothing
        E += Diagnostics.kinetic_energy_of_periodic_velocity_field(Ls) do x⃗
            external_fields.velocity(x⃗, t)
        end
    end
    E
end

function Diagnostics.kinetic_energy_nonperiodic(iter::VortexFilamentSolver; kws...)
    (; vs, fs,) = iter  # note: here we want vs (self-induced velocity) and not vL
    Ls = BiotSavart.periods(iter.prob.p)
    BiotSavart.domain_is_periodic(iter.prob.p) &&
        @warn(lazy"`kinetic_energy_nonperiodic` should only be called when working with non-periodic domains (got Ls = $Ls)")
    Diagnostics.kinetic_energy_nonperiodic(fs, vs, iter.prob.p; kws...)
end

function Diagnostics.energy_injection_rate(iter::VortexFilamentSolver, vL = iter.vL; kws...)
    (; vs, fs,) = iter
    p = iter.prob.p
    Diagnostics.energy_injection_rate(fs, vL, vs, p; kws...)
end

function Diagnostics.energy_flux(iter::VortexFilamentSolver, Nk_or_ks; kws...)
    velocities = (
        vs = (field = iter.vs, sign = -1),
        vinf = (field = CurvatureVector(), sign = -1),
    )
    p = iter.prob.p
    if hasproperty(iter, :vf)
        velocities = (; velocities..., vf = (field = iter.vf, sign = +1))
    end
    if hasproperty(iter, :vdiss)
        velocities = (; velocities..., vdiss = (field = iter.vdiss, sign = -1))
    end
    vs_buf = similar(iter.vs)
    Diagnostics.energy_flux(iter, iter.fs, velocities, Nk_or_ks, p; vs_buf, kws...)
end

function Diagnostics.energy_transfer_matrix(iter::VortexFilamentSolver, Nk_or_ks; kws...)
    params = iter.prob.p
    Diagnostics.energy_transfer_matrix(
        iter, iter.fs, iter.vs, Nk_or_ks, params; kws...
    )
end

# Note: filament_length is actually defined in the Filaments module, but we extend its
# alias in Diagnostics just for consistency with the other diagnostics (it doesn't make any
# difference really!).
function Diagnostics.filament_length(iter::VortexFilamentSolver; kws...)
    Diagnostics.filament_length(iter.fs; kws...)
end

function Diagnostics.vortex_impulse(iter::VortexFilamentSolver; kws...)
    Diagnostics.vortex_impulse(iter.fs; kws...)
end

function Diagnostics.helicity(iter::VortexFilamentSolver; kws...)
    Diagnostics.helicity(iter.fs, iter.vL, iter.prob.p; kws...)
end

function Diagnostics.stretching_rate(iter::VortexFilamentSolver; kws...)
    (; fs, vL,) = iter
    Diagnostics.stretching_rate(fs, vL; kws...)
end

# This allows passing a VortexFilamentSolver to energy_spectrum / energy_spectrum!.
Diagnostics.get_long_range_cache(iter::VortexFilamentSolver) = iter.cache_bs.longrange

function Diagnostics.integral_lengthscale(
        ks::AbstractVector{T}, Ek::AbstractVector{T}, Etot::T, Lvort::T, iter::VortexFilamentSolver{T}
    ) where {T}
    Diagnostics.integral_lengthscale(ks, Ek, Etot, Lvort, iter.prob.p)
end
