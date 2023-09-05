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

    tsub = t
    Mfast = scheme.nsubsteps
    cdt = dt / s
    hfast = cdt / Mfast  # timestep for evolution of fast component

    # Compute slow component at beginning of step
    rhs!(vS[1], fs, t, iter; component = Val(:fast))
    @. vS[1] = vs - vS[1]  # slow component at stage 1
    copy!(ftmp, fs)        # initial condition for stage 1

    # Always advect from latest filament positions in `ftmp`.
    fbase = ftmp

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

    # Stage 1. Note that the "slow" velocity is constant throughout this stage (we don't
    # need to introduce the normalised time τ).
    let i = 1
        # Solve auxiliary ODE in [t, t + dt/s].
        @assert tsub ≈ ts[i]
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))  # compute fast component
            @. vs = vs + vS[1]  # add slow component
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase)
            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vS[i + 1], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 2
    let i = 2
        # Solve auxiliary ODE in [t + dt/s, t + 2dt/s].
        @assert tsub ≈ ts[i]
        vslow = ntuple(j -> vS[j], Val(i))
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))
            τ = (t - ts[i]) / cdt  # normalised time in [0, 1]
            _MRI_inner_add_forcing_term!(vs, vslow, Γs, τ)
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase)
            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vS[i + 1], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 3
    let i = 3
        # Solve auxiliary ODE in [t + (i - 1) * dt/s, t + i * dt/s].
        @assert tsub ≈ ts[i]
        vslow = ntuple(j -> vS[j], Val(i))
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))
            τ = (t - ts[i]) / cdt  # normalised time in [0, 1]
            _MRI_inner_add_forcing_term!(vs, vslow, Γs, τ)
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase)
            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vS[i + 1], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 4
    let i = 4
        # Solve auxiliary ODE in [t + (i - 1) * dt/s, t + i * dt/s].
        @assert tsub ≈ ts[i]
        vslow = ntuple(j -> vS[j], Val(i))
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))
            τ = (t - ts[i]) / cdt  # normalised time in [0, 1]
            _MRI_inner_add_forcing_term!(vs, vslow, Γs, τ)
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase)
            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vS[i + 1], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 5
    let i = 5
        # Solve auxiliary ODE in [t + (i - 1) * dt/s, t + i * dt/s].
        @assert tsub ≈ ts[i]
        vslow = ntuple(j -> vS[j], Val(i))
        function rhs_inner!(vs, fs, t, iter)
            rhs!(vs, fs, t, iter; component = Val(:fast))
            τ = (t - ts[i]) / cdt  # normalised time in [0, 1]
            _MRI_inner_add_forcing_term!(vs, vslow, Γs, τ)
        end
        for m ∈ 1:Mfast
            rhs_inner!(vs, ftmp, tsub, iter)
            update_velocities!(
                rhs_inner!, advect!, cache_inner, iter;
                resize_cache = false, t = tsub, dt = hfast, fs = ftmp, vs = vs,
            )
            advect!(ftmp, vs, hfast; fbase)
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
