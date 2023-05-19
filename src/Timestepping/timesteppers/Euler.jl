"""
    Euler <: ExplicitTemporalScheme

Standard first-order Euler scheme.
"""
struct Euler <: ExplicitTemporalScheme end

# Euler doesn't need a cache...
struct EulerCache <: TemporalSchemeCache end

scheme(::EulerCache) = Euler()

function init_cache(::Euler, fs::VectorOfFilaments, vs::VectorOfArray)
    EulerCache()
end

function _update_velocities!(
        rhs!::F, advect!::G, cache::EulerCache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs, dt, to,) = iter
    rhs!(vs, fs, iter)
    vs
end
