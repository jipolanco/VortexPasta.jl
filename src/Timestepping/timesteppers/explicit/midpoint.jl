export Midpoint

"""
    Midpoint <: ExplicitScheme

Second-order, two-stage explicit midpoint method.
"""
struct Midpoint <: ExplicitScheme end

nbuf_filaments(::Midpoint) = 1
nbuf_velocities(::Midpoint) = 0

function _update_velocities!(
        ::Midpoint, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]

    # Stage 2
    advect!(ftmp, vs, dt/2; fbase = fs)
    rhs!(vs, ftmp, t + dt/2, iter)

    vs
end
