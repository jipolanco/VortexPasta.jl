export RK4

"""
    RK4 <: ExplicitScheme

Classic 4-stage Runge–Kutta method.
"""
struct RK4 <: ExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::RK4) = 1
nbuf_velocities(::RK4) = 1

function _update_velocities!(
        ::RK4, vs, rhs!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = iter.fs,
    ) where {F <: Function, G <: Function}
    (; fc, vc,) = cache
    ftmp = fc[1]
    vtmp = vc[1]

    # We assume that `vs` already contains the velocity at stage 1 (i.e. at the
    # current timestep). In other words, if I do `rhs!(vs, fs, iter)`, then
    # `vs` will have the same values it had before calling `rhs!`.

    # Stage 2
    advect!(ftmp, vs, dt/2; fbase = fs)
    rhs!(vtmp, ftmp, t + dt/2, iter)
    @. vs = vs + 2 * vtmp

    # Stage 3
    advect!(ftmp, vtmp, dt/2; fbase = fs)
    rhs!(vtmp, ftmp, t + dt/2, iter)
    @. vs = vs + 2 * vtmp

    # Stage 4 
    advect!(ftmp, vtmp, dt; fbase = fs)
    rhs!(vtmp, ftmp, t + dt, iter)

    # Final advecting velocity: v = (v[1] + 2 * v[2] + 2 * v[3] + v[4]) / 6
    @. vs = (vs + vtmp) / 6

    vs
end
