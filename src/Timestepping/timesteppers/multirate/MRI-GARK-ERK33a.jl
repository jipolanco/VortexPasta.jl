export SanduMRI33a

"""
    SanduMRI33a(nsubsteps::Int) <: MultirateScheme

3rd order, 3-stage multirate infinitesimal, generalised additive Runge–Kutta (MRI-GARK) scheme.

Uses a second order midpoint method for the fast component, with `nsubsteps` fast steps for
each RK stage.

This is the MRI-GARK-ERK33a method from Sandu, SIAM J. Numer. Anal. 57 (2019).
"""
struct SanduMRI33a <: MultirateScheme
    nsubsteps :: Int
end

nbuf_filaments(::SanduMRI33a) = 2
nbuf_velocities(::SanduMRI33a) = 4

function _update_velocities!(
        scheme::SanduMRI33a, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    fmid = fc[2]   # for midpoint method used for the fast component
    vS = ntuple(j -> vc[j], Val(3))  # slow velocity at each stage
    vslow = vc[4]  # total slow velocity

    tsub = t
    Mfast = scheme.nsubsteps
    cdt = dt / 3
    hfast = cdt / Mfast  # timestep for evolution of fast component

    # Compute slow component at beginning of step
    rhs!(vS[1], fs, t, iter; component = Val(:fast))
    @. vS[1] = vs - vS[1]  # slow component at stage 1
    copy!(ftmp, fs)        # initial condition for stage 1

    fbase = ftmp

    # Coupling coefficients (divided by Δc = 1/3)
    δ = -0.5  # completely empirical value (seems to work better than δ = -1/2 proposed by Sandu)
    Γ₀ = 3 * SMatrix{3, 3}(
        1/3, (-6δ - 7)/12, 0,          # column 1
        0, (6δ + 11)/12, (6δ - 5)/12,  # column 2
        0, 0, (3 - 2δ)/4,              # column 3
    )
    Γ₁ = 3 * SMatrix{3, 3}(
        0, (2δ + 1)/2, 1/2,
        0, -(2δ + 1)/2, 0,
        0, 0, δ,
    )

    ts = (t, t + cdt, t + 2cdt, t + dt)

    # Stage 1
    let i = 1, vslow = vS[1], fbase = ftmp
        # Solve auxiliary ODE in [t, t + dt/3].
        @assert tsub ≈ ts[i]
        for m ∈ 1:Mfast
            # Midpoint stage 1/2
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            @. vs = vs + vslow  # total velocity at this stage
            advect!(fmid, vs, hfast/2; fbase)

            # Midpoint stage 2/2
            rhs!(vs, fmid, tsub + hfast/2, iter; component = Val(:fast))
            @. vs = vs + vslow
            advect!(ftmp, vs, hfast; fbase)


            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vS[i + 1], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 2
    let i = 2, fbase = ftmp
        # Solve auxiliary ODE in [t + dt/3, t + 2dt/3].
        @assert tsub ≈ ts[i]
        for m ∈ 1:Mfast
            # Midpoint stage 1/2
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            τ = (tsub - ts[i]) / cdt  # normalised time in [0, 1]
            @. vs = vs + (
                (Γ₀[i, 1] + τ * Γ₁[i, 1]) * vS[1]
              + (Γ₀[i, 2] + τ * Γ₁[i, 2]) * vS[2]
            )
            advect!(fmid, vs, hfast/2; fbase)

            # Midpoint stage 2/2
            rhs!(vs, fmid, tsub + hfast/2, iter; component = Val(:fast))
            τ = (tsub + hfast/2 - ts[i]) / cdt
            @. vs = vs + (
                (Γ₀[i, 1] + τ * Γ₁[i, 1]) * vS[1]
              + (Γ₀[i, 2] + τ * Γ₁[i, 2]) * vS[2]
            )
            advect!(ftmp, vs, hfast; fbase)

            tsub += hfast
        end

        # Compute slow velocity at next stage
        rhs!(vS[i + 1], ftmp, tsub, iter; component = Val(:slow))
    end

    # Stage 3
    let i = 3
        # Solve auxiliary ODE in [t + 2dt/3, t + dt].
        @assert tsub ≈ ts[i]
        for m ∈ 1:Mfast
            # Midpoint stage 1/2
            rhs!(vs, ftmp, tsub, iter; component = Val(:fast))  # fast velocity at beginning of substep
            τ = (tsub - ts[i]) / cdt  # normalised time in [0, 1]
            @. vs = vs + (  # total slow velocity
                (Γ₀[i, 1] + τ * Γ₁[i, 1]) * vS[1]
              + (Γ₀[i, 2] + τ * Γ₁[i, 2]) * vS[2]
              + (Γ₀[i, 3] + τ * Γ₁[i, 3]) * vS[3]
            )
            advect!(fmid, vs, hfast/2; fbase)

            # Midpoint stage 2/2
            rhs!(vs, fmid, tsub + hfast/2, iter; component = Val(:fast))
            τ = (tsub + hfast/2 - ts[i]) / cdt
            @. vs = vs + (
                (Γ₀[i, 1] + τ * Γ₁[i, 1]) * vS[1]
              + (Γ₀[i, 2] + τ * Γ₁[i, 2]) * vS[2]
              + (Γ₀[i, 3] + τ * Γ₁[i, 3]) * vS[3]
            )
            advect!(ftmp, vs, hfast; fbase)

            tsub += hfast
        end
    end

    @assert tsub ≈ ts[4]

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end

function _inner_midpoint!(advect!::F, rhs!::G, fs, fbuf, vs, vslow, t, h, iter) where {F, G}
    @. vs = vs + vslow  # total velocity at this stage
    advect!(fbuf, vs, h/2; fbase = fs)
    rhs!(vs, fbuf, t + h/2, iter; component = Val(:fast))
    @. vs = vs + vslow
    advect!(fs, vs, h; fbase = fs)
    fs
end
