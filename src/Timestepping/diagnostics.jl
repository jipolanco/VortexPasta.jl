# Overload diagnostics functions for convenience, so that we can simply pass the current
# state of the solver to get diagnostics.
using ..Diagnostics: Diagnostics

function Diagnostics.kinetic_energy_from_streamfunction(iter::VortexFilamentSolver; kws...)
    (; ψs, fs,) = iter
    (; Γ, Ls,) = iter.prob.p.common
    Diagnostics.kinetic_energy_from_streamfunction(ψs, fs, Γ, Ls; kws...)
end

function Diagnostics.filament_length(iter::VortexFilamentSolver; kws...)
    Diagnostics.filament_length(iter.fs; kws...)
end

function Diagnostics.vortex_impulse(iter::VortexFilamentSolver; kws...)
    Diagnostics.vortex_impulse(iter.fs; kws...)
end

function Diagnostics.energy_spectrum(iter::VortexFilamentSolver; kws...)
    Diagnostics.energy_spectrum(iter.cache_bs.longrange; kws...)
end

function Diagnostics.energy_spectrum!(
        Ek::AbstractVector, ks::AbstractVector, iter::VortexFilamentSolver; kws...,
    )
    cache = iter.cache_bs.longrange
    Diagnostics.energy_spectrum!(Ek, ks, cache; kws...)
end
