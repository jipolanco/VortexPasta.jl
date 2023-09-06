export KenCarp4

"""
    KenCarp4 <: ImplicitExplicitScheme

4th order, 6 stage IMEX Runge–Kutta scheme by Kennedy & Carpenter (Appl. Numer. Math., 2003).
"""
struct KenCarp4 <: ImplicitExplicitScheme end

nstages(::KenCarp4) = 6

nbuf_filaments(::KenCarp4) = 1
nbuf_velocities(scheme::KenCarp4) = 2 * nstages(scheme)

function _update_velocities!(
        scheme::KenCarp4, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)
    s = nstages(scheme)

    ftmp = fc[1]
    vI = ntuple(j -> vc[0 + j], Val(s))      # implicit velocity at each stage
    vE = ntuple(j -> vc[s + j], Val(s - 1))  # explicit velocity at each stage
    v_explicit = vc[2s]

    tabE = tableau_explicit(scheme)
    tabI = tableau_implicit(scheme)
    (; bs, cs,) = tabE  # note: these vectors are identical in the implicit tableau
    AE = tabE.A
    AI = tabI.A

    # Stage 1
    let i = 1
        rhs!(vI[i], fs, t, iter; component = Val(:fast))  # compute fast component at stage 1
        @. vE[i] = vs - vI[i]  # slow component at stage 1
    end

    # Stage 2
    let i = 2, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase)  # initial guess for locations at stage i
        @. v_explicit = AE[i, 1] * vE[1] + AI[i, 1] * vI[1]  # known part of the velocity
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Separately compute both components at the final location
        rhs!(vI[i], ftmp, t + cdt, iter; component = Val(:fast))
        rhs!(vE[i], ftmp, t + cdt, iter; component = Val(:slow))
    end

    # Stage 3
    let i = 3, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage i
        @. v_explicit =  # known part of the velocity
            AE[i, 1] * vE[1] +
            AI[i, 1] * vI[1] +
            AE[i, 2] * vE[2] +
            AI[i, 2] * vI[2]
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Separately compute both components at the final location
        rhs!(vI[i], ftmp, t + cdt, iter; component = Val(:fast))
        rhs!(vE[i], ftmp, t + cdt, iter; component = Val(:slow))
    end

    # Stage 4
    let i = 4, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage i
        @. v_explicit =  # known part of the velocity
            AE[i, 1] * vE[1] +
            AI[i, 1] * vI[1] +
            AE[i, 2] * vE[2] +
            AI[i, 2] * vI[2] +
            AE[i, 3] * vE[3] +
            AI[i, 3] * vI[3]
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Separately compute both components at the final location
        rhs!(vI[i], ftmp, t + cdt, iter; component = Val(:fast))
        rhs!(vE[i], ftmp, t + cdt, iter; component = Val(:slow))
    end

    # Stage 5
    let i = 5, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage i
        @. v_explicit =  # known part of the velocity
            AE[i, 1] * vE[1] +
            AI[i, 1] * vI[1] +
            AE[i, 2] * vE[2] +
            AI[i, 2] * vI[2] +
            AE[i, 3] * vE[3] +
            AI[i, 3] * vI[3] +
            AE[i, 4] * vE[4] +
            AI[i, 4] * vI[4]
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Separately compute both components at the final location
        rhs!(vI[i], ftmp, t + cdt, iter; component = Val(:fast))
        rhs!(vE[i], ftmp, t + cdt, iter; component = Val(:slow))
    end

    # Stage 6
    let i = 6, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage i
        @. v_explicit =  # known part of the velocity
            AE[i, 1] * vE[1] +
            AI[i, 1] * vI[1] +
            AE[i, 2] * vE[2] +
            AI[i, 2] * vI[2] +
            AE[i, 3] * vE[3] +
            AI[i, 3] * vI[3] +
            AE[i, 4] * vE[4] +
            AI[i, 4] * vI[4] +
            AE[i, 5] * vE[5] +
            AI[i, 5] * vI[5]
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Since this is the final stage, we compute the full velocity at this stage onto vs
        # (we don't need the separate components)
        rhs!(vs, ftmp, t + cdt, iter; component = Val(:full))
    end

    # Final advecting velocity
    @. vs =
        bs[1] * (vE[1] + vI[1]) +
        bs[2] * (vE[2] + vI[2]) +
        bs[3] * (vE[3] + vI[3]) +
        bs[4] * (vE[4] + vI[4]) +
        bs[5] * (vE[5] + vI[5]) +
        bs[6] * vs  # vs actually contains the full velocity at stage 6

    vs
end

# See Appendix C of KC2003.
function tableau_explicit(::KenCarp4)
    cs = SVector(
        0,
        1/2,
        83/250,
        31/50,
        17/20,
        1,
    )
    bs = SVector(
        82889/524892,
        0,
        15625/83664,
        69875/102672,
        -2260/8211,
        1/4,
    )
    A = SMatrix{6, 6}(
        # Column 1
        0,
        1/2,
        13861/62500,
        -116923316275/2393684061468,
        -451086348788/2902428689909,
        647845179188/3216320057751,
        # Column 2
        0, 0,
        6889/62500,
        -2731218467317/15368042101831,
        -2682348792572/7519795681897,
        73281519250/8382639484533,
        # Column 3
        0, 0, 0,
        9408046702089/11113171139209,
        12662868775082/11960479115383,
        552539513391/3454668386233,
        # Column 4
        0, 0, 0, 0,
        3355817975965/11060851509271,
        3354512671639/8306763924573,
        # Column 5
        0, 0, 0, 0, 0,
        4040/17871,
        # Column 6
        0, 0, 0, 0, 0, 0,
    )
    (; A, bs, cs,)
end

# Note: this is an ESDIRK scheme (explicit, singly diagonally implicit Runge-Kutta)
function tableau_implicit(::KenCarp4)
    cs = SVector(
        0,
        1/2,
        83/250,
        31/50,
        17/20,
        1,
    )
    bs = SVector(
        82889/524892,
        0,
        15625/83664,
        69875/102672,
        -2260/8211,
        1/4,
    )
    A = SMatrix{6, 6}(
        # Column 1
        0,
        1/4,
        8611/62500,
        5012029/34652500,
        15267082809/155376265600,
        82889/524892,
        # Column 2
        0,
        1/4,
        -1743/31250,
        -654441/2922500,
        -71443401/120774400,
        0,
        # Column 3
        0, 0,
        1/4,
        174375/388108,
        730878875/902184768,
        15625/83664,
        # Column 4
        0, 0, 0,
        1/4,
        2285395/8070912,
        69875/102672,
        # Column 5
        0, 0, 0, 0,
        1/4,
        -2260/8211,
        # Column 6
        0, 0, 0, 0, 0,
        1/4,
    )
    (; A, bs, cs,)
end
