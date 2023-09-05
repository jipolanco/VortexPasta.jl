export SanduMRI45a

"""
    SanduMRI45a([inner = SSPRK33()], M::Int) <: MultirateScheme

4th order, 5-stage multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

Uses the explicit RK scheme `inner` for the fast component, with `M` inner steps for each
outer RK stage.

This is the MRI-GARK-ERK45 method from Sandu, SIAM J. Numer. Anal. 57 (2019).
"""
struct SanduMRI45a{InnerScheme <: ExplicitScheme} <: MultirateScheme
    inner     :: InnerScheme
    nsubsteps :: Int
end

SanduMRI45a(M::Integer) = SanduMRI45a(SSPRK33(), M)

nstages(::SanduMRI45a) = 5

nbuf_filaments(scheme::SanduMRI45a) = 1 + nbuf_filaments(inner_scheme(scheme))
nbuf_velocities(scheme::SanduMRI45a) = nstages(scheme) + nbuf_velocities(inner_scheme(scheme))

function _update_velocities!(
        scheme::SanduMRI45a, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    s = nstages(scheme)
    ftmp = fc[1]
    vS = ntuple(j -> vc[j], Val(s))  # slow velocity at each stage

    cache_inner = let scheme = inner_scheme(scheme)
        TemporalSchemeCache(
            scheme,
            ntuple(j -> fc[1 + j], Val(nbuf_filaments(scheme))),
            ntuple(j -> vc[s + j], Val(nbuf_velocities(scheme))),
        )
    end

    Mfast = scheme.nsubsteps
    cdt = dt / s
    hfast = cdt / Mfast  # timestep for evolution of fast component

    # Compute slow component at beginning of step
    rhs!(vS[1], fs, t, iter; component = Val(:fast))
    @. vS[1] = vs - vS[1]  # slow component at stage 1
    copy!(ftmp, fs)        # initial condition for stage 1

    # Coupling coefficients (divided by Δc = 1/5)
    Γ₀ = 5 * SMatrix{5, 5}(
        1/5, -53/16, -36562993/71394880, -7631593/71394880, 277061/303808,  # column 1
        0, 281/80, 34903117/17848720, -166232021/35697440, -209323/1139280,
        0, 0, -88770499/71394880, 6068517/1519040, -1360217/1139280,
        0, 0, 0, 8644289/8924360, -148789/56964,
        0, 0, 0, 0, 147889/45120,
    )
    Γ₁ = 5 * SMatrix{5, 5}(
        0, 503/80, -1365537/35697440, 66974357/35697440, -18227/7520,
        0, -503/80, 4963773/7139488, 21445367/7139488, 2,
        0, 0, -1465833/2231090, -3, 1,
        0, 0, 0, -8388609/4462180, 5,
        0, 0, 0, 0, -41933/7520,
    )
    Γs = (Γ₀, Γ₁)

    ts = ntuple(j -> t + (j - 1) * cdt, Val(s + 1))
    @assert ts[s + 1] ≈ t + dt

    let i = 1
        _MRI_stage!(Val(i), rhs!, advect!, ftmp, vS, ts, iter, vs, cache_inner, Mfast, hfast, Γs)
    end
    let i = 2
        rhs!(vS[i], ftmp, ts[i], iter; component = Val(:slow))
        _MRI_stage!(Val(i), rhs!, advect!, ftmp, vS, ts, iter, vs, cache_inner, Mfast, hfast, Γs)
    end
    let i = 3
        rhs!(vS[i], ftmp, ts[i], iter; component = Val(:slow))
        _MRI_stage!(Val(i), rhs!, advect!, ftmp, vS, ts, iter, vs, cache_inner, Mfast, hfast, Γs)
    end
    let i = 4
        rhs!(vS[i], ftmp, ts[i], iter; component = Val(:slow))
        _MRI_stage!(Val(i), rhs!, advect!, ftmp, vS, ts, iter, vs, cache_inner, Mfast, hfast, Γs)
    end
    let i = 5
        rhs!(vS[i], ftmp, ts[i], iter; component = Val(:slow))
        _MRI_stage!(Val(i), rhs!, advect!, ftmp, vS, ts, iter, vs, cache_inner, Mfast, hfast, Γs)
    end

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end
