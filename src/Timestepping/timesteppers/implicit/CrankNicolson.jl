export CrankNicolson

"""
    CrankNicolson() <: ImplicitScheme

2nd order Crank–Nicolson implicit scheme.
"""
struct CrankNicolson <: ImplicitScheme end

nbuf_filaments(::CrankNicolson) = 1
nbuf_velocities(::CrankNicolson) = 1

function _update_velocities!(
        ::CrankNicolson, vs, rhs!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = iter.fs,
    ) where {F <: Function, G <: Function}
    (; fc, vc,) = cache
    ftmp = fc[1]
    v_explicit = vc[1]

    @. v_explicit = vs / 2  # explicit part of the RHS

    # Initial guess for locations at time t + dt
    #
    # - Explicit Euler
    # advect!(ftmp, vs, dt; fbase = fs)
    #
    # - Explicit midpoint (RK2) -- helps with convergence!
    advect!(ftmp, vs, dt/2; fbase = fs)
    rhs!(vs, ftmp, t + dt/2, iter)
    advect!(ftmp, vs, dt; fbase = fs)

    # At the end, `vs` contains the advecting velocity = (v[n] + v[n + 1]) / 2.
    solve_fixed_point!(
        ftmp, rhs!, advect!, iter, vs, v_explicit;
        t, dt, fbase = fs, aI_diag = 1/2,
    )

    vs
end
