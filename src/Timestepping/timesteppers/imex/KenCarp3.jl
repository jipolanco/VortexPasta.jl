export KenCarp3

"""
    KenCarp3 <: ImplicitExplicitScheme

3rd order, 4 stage IMEX Runge–Kutta scheme by Kennedy & Carpenter (Appl. Numer. Math., 2003).
"""
struct KenCarp3 <: ImplicitExplicitScheme end

nstages(::KenCarp3) = 4

nbuf_filaments(::KenCarp3) = 1
nbuf_velocities(scheme::KenCarp3) = 2 * nstages(scheme)

function _update_velocities!(
        scheme::KenCarp3, rhs!::F, advect!::G, cache, iter::AbstractSolver,
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

    # Implicit part of the RHS.
    rhs_implicit! = imex_rhs_implicit(rhs!)
    rhs_explicit! = imex_rhs_explicit(rhs!)
    rhs_full! = imex_rhs_full(rhs!)

    # Stage 1
    let i = 1
        rhs_implicit!(vI[i], fs, t, iter)  # compute fast component at stage 1
        @. vE[i] = vs - vI[i]  # slow component at stage 1
    end

    # Stage 2
    let i = 2, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase)  # initial guess for locations at stage i
        _imex_explicit_rhs!(Val(i), v_explicit, (AE, AI), (vE, vI))  # known part of the velocity
        solve_fixed_point!(
            ftmp, rhs_implicit!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Separately compute both components at the final location
        rhs_implicit!(vI[i], ftmp, t + cdt, iter)
        rhs_explicit!(vE[i], ftmp, t + cdt, iter)
    end

    # Stage 3
    let i = 3, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage i
        _imex_explicit_rhs!(Val(i), v_explicit, (AE, AI), (vE, vI))  # known part of the velocity
        solve_fixed_point!(
            ftmp, rhs_implicit!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Separately compute both components at the final location
        rhs_implicit!(vI[i], ftmp, t + cdt, iter)
        rhs_explicit!(vE[i], ftmp, t + cdt, iter)
    end

    # Stage 4
    let i = 4, fbase = fs, vtmp = vI[i]
        cdt = cs[i] * dt
        advect!(ftmp, vs, cdt; fbase = fs)  # initial guess for locations at stage i
        _imex_explicit_rhs!(Val(i), v_explicit, (AE, AI), (vE, vI))  # known part of the velocity
        solve_fixed_point!(
            ftmp, rhs_implicit!, advect!, iter, vtmp, v_explicit;
            cdt, fbase, aI_diag = AI[i, i],
        )
        # Since this is the final stage, we compute the full velocity at this stage onto vs
        # (we don't need the separate components)
        rhs_full!(vs, ftmp, t + cdt, iter)
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
