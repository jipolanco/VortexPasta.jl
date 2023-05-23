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
    (; vs,) = iter
    # We assume that `vs` already contains the velocity at the current
    # timestep. In other words, if I do `rhs!(vs, fs, iter)`, then `vs` will
    # have the same values it had before calling `rhs!`.
    vs
end
