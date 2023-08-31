export KenCarp3

"""
    KenCarp3 <: ImplicitExplicitScheme

3rd order IMEX Runge–Kutta scheme by Kennedy & Carpenter (Appl. Numer. Math., 2003).
"""
struct KenCarp3 <: ImplicitExplicitScheme end

nbuf_filaments(::KenCarp3) = 1
nbuf_velocities(::KenCarp3) = 8  # TODO reduce?

function _update_velocities!(
        scheme::KenCarp3, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    ftmp = fc[1]
    vE = ntuple(i -> vc[0 + i], Val(3))  # explicit velocity at each stage
    vI = ntuple(i -> vc[3 + i], Val(4))  # implicit velocity at each stage
    v_explicit = vc[8]

    tabE = tableau_explicit(scheme)
    tabI = tableau_implicit(scheme)
    (; bs, cs,) = tabE  # note: these vectors are identical in the implicit tableau
    AE = tabE.A
    AI = tabI.A

    # Parameters for fixed-point iterations.
    rtol = 1e-12
    nmax = 100

    # Stage 1
    let s = 1
        rhs!(vI[s], fs, t, iter; component = Val(:fast))  # compute fast component at stage 1
        @. vE[s] = vs - vI[s]  # slow component at stage 1 (note: vE is aliased to vs)
    end

    # Stage 2
    let s = 2, fbase = fs, vtmp = vI[s]
        cdt = cs[s] * dt
        advect!(ftmp, vs, cdt; fbase)  # initial guess for locations at stage s
        @. v_explicit = AE[s, 1] * vE[1] + AI[s, 1] * vI[1]  # known part of the velocity
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_ss = AI[s, s],
        )
        # Separately compute both components at the final location
        rhs!(vI[s], ftmp, t + cdt, iter; component = Val(:fast))
        rhs!(vE[s], ftmp, t + cdt, iter; component = Val(:slow))
    end

    # Stage 3
    let s = 3, fbase = fs, vtmp = vI[s]
        cdt = cs[s] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage s
        @. v_explicit =  # known part of the velocity
            AE[s, 1] * vE[1] +
            AE[s, 2] * vE[2] +
            AI[s, 1] * vI[1] +
            AI[s, 2] * vI[2]
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_ss = AI[s, s],
        )
        # Separately compute both components at the final location
        rhs!(vI[s], ftmp, t + cdt, iter; component = Val(:fast))
        rhs!(vE[s], ftmp, t + cdt, iter; component = Val(:slow))
    end

    # Stage 4
    let s = 4, fbase = fs, vtmp = vI[s]
        cdt = cs[s] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage s
        @. v_explicit =  # known part of the velocity
            AE[s, 1] * vE[1] +
            AE[s, 2] * vE[2] +
            AE[s, 3] * vE[3] +
            AI[s, 1] * vI[1] +
            AI[s, 2] * vI[2] +
            AI[s, 3] * vI[3]
        solve_fixed_point!(
            ftmp, rhs!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_ss = AI[s, s],
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
        bs[4] * vs  # vs actually contains the full velocity at stage 4

    vs
end

# See Appendix C of KC2003.
function tableau_explicit(::KenCarp3)
    cs = SVector(
        0.0,
        1767732205903/2027836641118,
        3/5,
        1.0,
    )
    bs = SVector(
        1471266399579/7840856788654,
        -4482444167858/7529755066697,
        11266239266428/11593286722821,
        1767732205903/4055673282236,
    )
    A = SMatrix{4, 4}(
        # Column 1
        0,
        1767732205903/2027836641118,
        5535828885825/10492691773637,
        6485989280629/16251701735622,
        # Column 2
        0, 0,
        788022342437/10882634858940,
        -4246266847089/9704473918619,
        # Column 3
        0, 0, 0,
        10755448449292/10357097424841,
        # Column 4
        0, 0, 0, 0,
    )
    (; A, bs, cs,)
end

# Note: this is an ESDIRK scheme (explicit, singly diagonally implicit Runge-Kutta)
function tableau_implicit(::KenCarp3)
    cs = SVector(
        0.0,
        1767732205903/2027836641118,
        3/5,
        1.0,
    )
    bs = SVector(
        1471266399579/7840856788654,
        -4482444167858/7529755066697,
        11266239266428/11593286722821,
        1767732205903/4055673282236,
    )
    A = SMatrix{4, 4}(
        # Column 1
        0,
        1767732205903/4055673282236,
        2746238789719/10658868560708,
        1471266399579/7840856788654,
        # Column 2
        0,
        1767732205903/4055673282236,
        -640167445237/6845629431997,
        -4482444167858/7529755066697,
        # Column 3
        0, 0,
        1767732205903/4055673282236,
        11266239266428/11593286722821,
        # Column 4
        0, 0, 0,
        1767732205903/4055673282236,
    )
    (; A, bs, cs,)
end
