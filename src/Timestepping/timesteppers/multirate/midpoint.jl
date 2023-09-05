export MultirateMidpoint

"""
    MultirateMidpoint([inner = Euler()], M::Int) <: MultirateScheme

2nd order, two-stage multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

Uses the explicit RK scheme `inner` for the fast component (by default a 1st-order Euler
scheme), with `M` inner steps for each outer RK stage.

This is the MRI-GARK-ERK22a method from Sandu, SIAM J. Numer. Anal. 57 (2019).
"""
struct MultirateMidpoint{InnerScheme <: ExplicitScheme} <: MultirateScheme
    inner     :: InnerScheme
    nsubsteps :: Int
end

MultirateMidpoint(M::Integer) = MultirateMidpoint(Euler(), M)

inner_scheme(scheme::MultirateMidpoint) = scheme.inner

nbuf_filaments(scheme::MultirateMidpoint) = 1 + nbuf_filaments(inner_scheme(scheme))
nbuf_velocities(scheme::MultirateMidpoint) = 1 + nbuf_velocities(inner_scheme(scheme))

function _update_velocities!(
        scheme::MultirateMidpoint, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    vslow = vc[1]  # "slow" velocity at each stage

    cache_inner = let scheme = inner_scheme(scheme)
        TemporalSchemeCache(
            scheme,
            ntuple(j -> fc[1 + j], Val(nbuf_filaments(scheme))),
            ntuple(j -> vc[1 + j], Val(nbuf_velocities(scheme))),
        )
    end

    tsub = t    # current time of ftmp
    Mfast = scheme.nsubsteps  # number of Euler substeps
    hfast = (dt/2) / Mfast  # timestep for evolution of fast component

    # Compute slow component at beginning of step
    rhs!(vslow, fs, t, iter; component = Val(:fast))
    @. vslow = vs - vslow  # slow component at stage 1
    copy!(ftmp, fs)        # initial condition for stage 1

    # Stage 1
    let i = 1
        # Solve auxiliary ODE in [t, t + dt/2].
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))
            @. vs = vs + vslow
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase = ftmp)
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
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))
            @. vs = vs + vslow
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase = ftmp)
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
