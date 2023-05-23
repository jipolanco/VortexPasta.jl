"""
    SSPRK33 <: ExplicitTemporalScheme

Three-stage, third-order strong stability preserving (SSP) method.

See <https://en.wikipedia.org/wiki/List_of_Runge%E2%80%93Kutta_methods#Third-order_Strong_Stability_Preserving_Runge-Kutta_(SSPRK3)>.
"""
struct SSPRK33 <: ExplicitTemporalScheme end

struct SSPRK33Cache{
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfArray{<:Vec3},
    } <: TemporalSchemeCache
    fc :: Tuple{Filaments}
    vc :: Tuple{Velocities}
end

scheme(::SSPRK33Cache) = SSPRK33()

function init_cache(::SSPRK33, fs::VectorOfFilaments, vs::VectorOfArray)
    fc = (map(similar, fs),)
    vc = (similar(vs),)
    SSPRK33Cache(fc, vc)
end

function _update_velocities!(
        rhs!::F, advect!::G, cache::SSPRK33Cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs, t, dt, to,) = iter
    (; fc, vc,) = cache

    resize!(cache, fs)  # in case the number of nodes (or filaments) has changed
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
