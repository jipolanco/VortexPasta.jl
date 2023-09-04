export MultirateMidpoint

"""
    MultirateMidpoint(nsubsteps::Int = 16) <: MultirateScheme

2nd order multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

This is the MRI-GARK-ERK22a method from Sandu, SIAM J. Numer. Anal. 57 (2019).
"""
struct MultirateMidpoint <: MultirateScheme
    nsubsteps :: Int
end

nbuf_filaments(::MultirateMidpoint) = 1
nbuf_velocities(::MultirateMidpoint) = 2

function _update_velocities!(
        scheme::MultirateMidpoint, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    vS = ntuple(j -> vc[j], Val(2))  # "slow" velocity at each stage
    vF = vs  # reuse memory for "fast" velocity

    tsub = t    # current time of ftmp
    M = scheme.nsubsteps  # number of Euler substeps

    # Stage 1
    let i = 1
        # Solve auxiliary ODE in [t, t + dt/2].
        # Use Euler scheme with very small timestep.
        h = (dt/2) / M

        # Initial advection: we start from `fs` and use the total velocity at the
        # beginning of the timestep.
        rhs!(vS[2], fs, t, iter; component = Val(:fast))  # compute fast component at beginning of step (reusing vS[2])
        @. vS[1] = vF - vS[2]  # slow component at stage 1 (note: vF === vs)
        advect!(ftmp, vF, h; fbase = fs)
        tsub += h

        # Next advections: update "fast" velocity and advect from `ftmp`.
        for m ∈ 2:M
            rhs!(vF, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            @. vF = vF + vS[1] # total velocity at this step
            advect!(ftmp, vF, h; fbase = ftmp)
            tsub += h
        end
        @assert tsub ≈ t + dt/2

        # Now vF contains the fast component at the final time of the auxiliary ODE.
        # We want the slow component at that time.
        rhs!(vS[2], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 2
    let i = 2
        # Solve auxiliary ODE in [t + dt/2, t + dt].
        @assert tsub ≈ t + dt/2
        h = (dt/2) / M
        @. vS[2] = 2 * vS[2] - vS[1]  # total slow velocity in this range
        for m ∈ 1:M
            rhs!(vF, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            @. vF = vF + vS[2]  # fast + slow velocity
            advect!(ftmp, vF, h; fbase = ftmp)
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
