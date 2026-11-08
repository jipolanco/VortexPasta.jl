function Diagnostics.kinetic_energy(iter::PointVortexSolver)
    Diagnostics.kinetic_energy(iter.vortices, iter.quantities.ψs, get_params(iter))
end
