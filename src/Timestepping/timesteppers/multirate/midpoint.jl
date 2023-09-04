export MultirateMidpoint

"""
    MultirateMidpoint(nsubsteps::Int) <: MultirateScheme

2nd order, two-stage multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

Uses a first order Euler method for the fast component, with `nsubsteps` Euler substeps for
each RK stage.

This is the MRI-GARK-ERK22a method from Sandu, SIAM J. Numer. Anal. 57 (2019).
"""
struct MultirateMidpoint <: MultirateScheme
    nsubsteps :: Int
end

nbuf_filaments(::MultirateMidpoint) = 1
nbuf_velocities(::MultirateMidpoint) = 1

function _update_velocities!(
        scheme::MultirateMidpoint, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    vslow = vc[1]  # "slow" velocity at each stage

    tsub = t    # current time of ftmp
    M = scheme.nsubsteps  # number of Euler substeps

    # Stage 1
    let i = 1
        # Solve auxiliary ODE in [t, t + dt/2].
        # Use Euler scheme with very small timestep.
        h = (dt/2) / M

        # Initial advection: we start from `fs` and use the total velocity at the
        # beginning of the timestep.
        rhs!(vslow, fs, t, iter; component = Val(:fast))  # compute fast component at beginning of step (reusing vslow[2])
        @. vslow = vs - vslow  # slow component at stage 1
        let fbase = fs
            advect!(ftmp, vs, h; fbase)
        end
        tsub += h

        # Next advections: update "fast" velocity and advect from `ftmp`.
        for m ∈ 2:M
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            let fbase = ftmp
                @. vs = vs + vslow  # total velocity at this step
                advect!(ftmp, vs, h; fbase)
            end
            tsub += h
        end
        @assert tsub ≈ t + dt/2

        # Now vs contains the fast component at the final time of the auxiliary ODE.
        # We want the slow component at that time.
        rhs!(vs, ftmp, tsub, iter; component = Val(:slow))
        @. vslow = 2 * vs - vslow  # total slow velocity for stage 2
    end

    # Stage 2
    let i = 2
        # Solve auxiliary ODE in [t + dt/2, t + dt].
        @assert tsub ≈ t + dt/2
        h = (dt/2) / M
        for m ∈ 1:M
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            let fbase = ftmp
                @. vs = vs + vslow  # total velocity at this step
                advect!(ftmp, vs, h; fbase)
            end
            tsub += h
        end
        @assert tsub ≈ t + dt
    end

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end
