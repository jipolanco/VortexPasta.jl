export Euler

"""
    Euler <: ExplicitScheme

Standard first-order Euler scheme.
"""
struct Euler <: ExplicitScheme end

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::Euler) = 0
nbuf_velocities(::Euler) = 0

function _update_velocities!(
        ::Euler, rhs!::F, advect!::G, cache, iter::AbstractSolver;
        vs = iter.vs, kws...,
    ) where {F <: Function, G <: Function}
    # We assume that `vs` already contains the velocity at the current
    # timestep. In other words, if I do `rhs!(vs, fs, iter)`, then `vs` will
    # have the same values it had before calling `rhs!`.
    vs
end
