"""
    RK4 <: ExplicitTemporalScheme

Classic 4-step Runge–Kutta method.
"""
struct RK4 <: ExplicitTemporalScheme end

function init_cache(::RK4, fs::VectorOfFilaments, vs::VectorOfArray)
    fc = map(similar, fs) :: VectorOfFilaments
    vc = (similar(vs),)
    RK4Cache(fc, vc)
end

struct RK4Cache{
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfArray{<:Vec3},
    } <: TemporalSchemeCache
    fc :: Filaments
    vc :: NTuple{1, Velocities}
end

scheme(::RK4Cache) = RK4()

function _update_velocities!(
        rhs!::F, advect!::G, cache::RK4Cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs, dt, to,) = iter
    (; fc, vc,) = cache

    vtmp = vc[1]

    @assert length(vs) == length(fs)
    resize!(cache, fs)  # in case the number of nodes (or filaments) has changed

    # Step 1
    rhs!(vs, fs, iter)

    # Step 2
    advect!(fc, vs, dt/2; fbase = fs)
    rhs!(vtmp, fc, iter)
    @. vs = vs + 2 * vtmp

    # Step 3
    advect!(fc, vtmp, dt/2; fbase = fs)
    rhs!(vtmp, fc, iter)
    @. vs = vs + 2 * vtmp

    # Step 4 
    advect!(fc, vtmp, dt; fbase = fs)
    rhs!(vtmp, fc, iter)

    # Final velocity: v = (v[1] + 2 * v[2] + 2 * v[3] + v[4]) / 6
    @. vs = (vs + vtmp) / 6

    vs
end
