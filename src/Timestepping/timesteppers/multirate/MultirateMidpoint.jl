export MultirateMidpoint

"""
    MultirateMidpoint() <: MultirateScheme

2nd order multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

See Sandu, SIAM J. Numer. Anal. 57 (2019).
"""
struct MultirateMidpoint <: MultirateScheme end

nbuf_filaments(::MultirateMidpoint) = 1
nbuf_velocities(::MultirateMidpoint) = 3

function _update_velocities!(
        ::MultirateMidpoint, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    vS = ntuple(j -> vc[j], Val(2))  # slow velocity at each stage
    vF = vc[3]

    copy!(ftmp, fs)  # initial locations
    tsub = t  # current time of ftmp

    M = 32  # number of Euler substeps

    # Stage 1
    let i = 1
        rhs!(vF, fs, t, iter; component = Val(:fast))  # compute fast component at beginning of step
        @. vS[1] = vs - vF  # slow component at stage 1

        # Solve auxiliary ODE in [t, t + dt/2].
        # Use Euler scheme with very small timestep.
        h = (dt/2) / M
        for m ∈ 1:M
            @. vF = vF + vS[1] # total velocity at this step
            tsub += h
            advect!(ftmp, vF, h; fbase = ftmp)
            rhs!(vF, ftmp, tsub, iter; component = Val(:fast))
        end
        @assert tsub ≈ t + dt/2

        # Now vF contains the fast component at the final time of the auxiliary ODE.
        # We want the slow component at that time.
        rhs!(vS[2], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 2
    let i = 2, v_slow = vS[2]
        # Solve auxiliary ODE in [t + dt/2, t + dt].
        @assert tsub ≈ t + dt/2
        h = (dt/2) / M
        @. v_slow = 2 * vS[2] - vS[1]
        for m ∈ 1:M
            @. vF = vF + v_slow
            tsub += h
            advect!(ftmp, vF, h; fbase = ftmp)
            rhs!(vF, ftmp, tsub, iter; component = Val(:fast))  # TODO last evaluation is not needed!
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
