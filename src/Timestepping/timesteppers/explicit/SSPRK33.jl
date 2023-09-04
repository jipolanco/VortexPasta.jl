export SSPRK33

"""
    SSPRK33 <: ExplicitScheme

Three-stage, third-order strong stability preserving (SSP) method.

See [Wikipedia](https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Third-order_Strong_Stability_Preserving_Runge-Kutta_(SSPRK3))
for details.
"""
struct SSPRK33 <: ExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::SSPRK33) = 1
nbuf_velocities(::SSPRK33) = 1

function _update_velocities!(
        ::SSPRK33, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    vtmp = vc[1]

    # Stage 2
    advect!(ftmp, vs, dt; fbase = fs)
    rhs!(vtmp, ftmp, t + dt, iter)
    @. vs = vs + vtmp  # = v₁ + v₂

    # Stage 3
    advect!(ftmp, vs, dt/4; fbase = fs)  # same as advecting with vs/2 and dt/2
    rhs!(vtmp, ftmp, t + dt/2, iter)

    # Final advecting velocity: v = (v[1] + 2 * v[2] + 4 * v[3]) / 6
    @. vs = (vs + 4 * vtmp) / 6

    vs
end
