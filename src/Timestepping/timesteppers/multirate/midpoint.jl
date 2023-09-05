export MultirateMidpoint

"""
    MultirateMidpoint(nsubsteps::Int) <: MultirateScheme

2nd order, two-stage multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

Uses a first order Euler method for the fast component, with `nsubsteps` inner steps for
each outer RK stage.

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
    Mfast = scheme.nsubsteps  # number of Euler substeps
    hfast = (dt/2) / Mfast  # timestep for evolution of fast component

    # Compute slow component at beginning of step
    rhs!(vslow, fs, t, iter; component = Val(:fast))
    @. vslow = vs - vslow  # slow component at stage 1
    copy!(ftmp, fs)        # initial condition for stage 1

    fbase = ftmp

    # Stage 1
    let i = 1
        # Solve auxiliary ODE in [t, t + dt/2].
        for m ∈ 1:Mfast
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            @. vs = vs + vslow  # total velocity at this step
            advect!(ftmp, vs, hfast; fbase)
            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vs, ftmp, tsub, iter; component = Val(:slow))
        @. vslow = 2 * vs - vslow
    end

    @assert tsub ≈ t + dt/2

    # Stage 2
    let i = 2
        # Solve auxiliary ODE in [t + dt/2, t + dt].
        @assert tsub ≈ t + dt/2
        for m ∈ 1:Mfast
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            let fbase = ftmp
                @. vs = vs + vslow  # total velocity at this step
                advect!(ftmp, vs, hfast; fbase)
            end
            tsub += hfast
        end
    end

    @assert tsub ≈ t + dt

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end
