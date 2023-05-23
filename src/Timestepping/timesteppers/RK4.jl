"""
    RK4 <: ExplicitTemporalScheme

Classic 4-step Runge–Kutta method.
"""
struct RK4 <: ExplicitTemporalScheme end

struct RK4Cache{
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfArray{<:Vec3},
    } <: TemporalSchemeCache
    fc :: Tuple{Filaments}
    vc :: Tuple{Velocities}
end

scheme(::RK4Cache) = RK4()

function init_cache(::RK4, fs::VectorOfFilaments, vs::VectorOfArray)
    fc = (map(similar, fs),)
    vc = (similar(vs),)
    RK4Cache(fc, vc)
end

function _update_velocities!(
        rhs!::F, advect!::G, cache::RK4Cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs, dt, to,) = iter
    (; fc, vc,) = cache

    ftmp = fc[1]
    vtmp = vc[1]

    @assert length(vs) == length(fs)
    resize!(cache, fs)  # in case the number of nodes (or filaments) has changed

    # We assume that `vs` already contains the velocity at step 1 (i.e. at the
    # current timestep). In other words, if I do `rhs!(vs, fs, iter)`, then
    # `vs` will have the same values it had before calling `rhs!`.

    # Step 2
    t = iter.t + dt/2
    advect!(ftmp, vs, dt/2; fbase = fs)
    rhs!(vtmp, ftmp, t, iter)
    @. vs = vs + 2 * vtmp

    # Step 3
    t = iter.t + dt/2
    advect!(ftmp, vtmp, dt/2; fbase = fs)
    rhs!(vtmp, ftmp, t, iter)
    @. vs = vs + 2 * vtmp

    # Step 4 
    t = iter.t + dt
    advect!(ftmp, vtmp, dt; fbase = fs)
    rhs!(vtmp, ftmp, t, iter)

    # Final advecting velocity: v = (v[1] + 2 * v[2] + 2 * v[3] + v[4]) / 6
    @. vs = (vs + vtmp) / 6

    vs
end
