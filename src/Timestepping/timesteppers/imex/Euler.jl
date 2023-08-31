export IMEXEuler

"""
    IMEXEuler <: ImplicitExplicitScheme

Forward-backward Euler scheme.

This is the (1,2,1) scheme in the notation of Ascher et al. (1997).
"""
struct IMEXEuler <: ImplicitExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::IMEXEuler) = 1
nbuf_velocities(::IMEXEuler) = 1

function _update_velocities!(
        ::IMEXEuler, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]

    # Initial guess for locations at stage 2
    advect!(ftmp, vs, dt; fbase = fs)

    # We assume that `vs` contains the "full" velocity at stage 1 (i.e. at the current
    # timestep). We compute and subtract the fast component to obtain the slow component.
    vE = vs     # explicit (slow) component
    vI = vc[1]  # implicit (fast) component
    rhs!(vI, fs, t, iter; component = Val(:fast))
    @. vE = vs - vI  # slow component at stage 1 (note: vE is aliased to vs)

    solve_fixed_point!(
        ftmp, rhs!, advect!, iter, vI, vE;
        cdt = dt, fbase = fs, aI_diag = 1,
    )

    # Final advecting velocity
    t_final = t + dt
    rhs!(vs, ftmp, t_final, iter; component = Val(:full))

    vs
end
