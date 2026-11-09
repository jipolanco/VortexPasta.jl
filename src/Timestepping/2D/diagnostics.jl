function Diagnostics.kinetic_energy(iter::PointVortexSolver)
    Diagnostics.kinetic_energy(iter.vortices, iter.quantities.ψs, get_params(iter))
end

function Diagnostics.vortex_impulse(iter::PointVortexSolver)
    Diagnostics.vortex_impulse(iter.vortices)
end
