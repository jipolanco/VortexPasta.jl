function Diagnostics.kinetic_energy(iter::PointVortexSolver)
    _check_diagnostics(iter)
    Diagnostics.kinetic_energy(iter.vortices, iter.quantities.ψs, get_params(iter))
end

function Diagnostics.vortex_impulse(iter::PointVortexSolver)
    _check_diagnostics(iter)
    Diagnostics.vortex_impulse(iter.vortices)
end

function Diagnostics.energy_flux(iter::PointVortexSolver, Nk_or_ks)
    _check_diagnostics(iter)
    params = get_params(iter)
    velocities = (
        vs = (field = iter.vs, sign = -1),
    )
    if hasproperty(iter.quantities, :vf)
        velocities = (; velocities..., vf = (field = iter.vf, sign = +1))
    end
    if hasproperty(iter.quantities, :vdiss)
        velocities = (; velocities..., vdiss = (field = iter.vdiss, sign = -1))
    end
    vs_buf = similar(iter.vs)
    Diagnostics.energy_flux(iter, iter.vortices, velocities, Nk_or_ks, params; vs_buf)
end
