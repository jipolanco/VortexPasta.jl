# NOTE: this scheme has very similar properties to KenCarp3

export Ascher343

"""
    Ascher343 <: ImplicitExplicitScheme

3rd order, 4 stage IMEX Runge–Kutta scheme by Ascher et al. (Appl. Numer. Math., 1997).

⚠ This scheme may be removed in the future as it behaves very similarly to [`KenCarp3`](@ref).
"""
struct Ascher343 <: ImplicitExplicitScheme end

nstages(::Ascher343) = 4

nbuf_filaments(::Ascher343) = 1
nbuf_velocities(scheme::Ascher343) = 2 * nstages(scheme)

function _update_velocities!(
        scheme::Ascher343, rhs!::F, advect!::G, cache, iter::AbstractSolver,
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
        @. v_explicit = (  # known part of the velocity
            AE[i, 1] * vE[1]
          # + AI[i, 1] * vI[1]  # note: AI[:, 1] == 0
        )
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
        @. v_explicit = (  # known part of the velocity
            AE[i, 1] * vE[1]
          # + AI[i, 1] * vI[1]
          + AE[i, 2] * vE[2]
          + AI[i, 2] * vI[2]
        )
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
        @. v_explicit = (  # known part of the velocity
            AE[i, 1] * vE[1]
          # + AI[i, 1] * vI[1]
          + AE[i, 2] * vE[2]
          + AI[i, 2] * vI[2]
          + AE[i, 3] * vE[3]
          + AI[i, 3] * vI[3]
        )
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
        # FIXME inference fails when factorising the expression
        # Maybe it's a problem with broadcasting of VectorOfVector's?
        # bs[1] * (vE[1] + vI[1]) +
        bs[1] * vE[1] + bs[1] * vI[1] +
        bs[2] * vE[2] + bs[2] * vI[2] +
        bs[3] * vE[3] + bs[3] * vI[3] +
        bs[4] * vs  # vs actually contains the full velocity at stage 4

    vs
end

function tableau_explicit(::Ascher343)
    γ = 0.435866521508459
    b1 = -3/2 * γ^2 + 4γ - 1/4
    b2 = 1 - (γ + b1)
    cs = SVector(0.0, γ, (1 + γ) / 2, 1.0)
    bs = SVector(0, b1, b2, γ)
    a42 = 0.5529291479
    a43 = a42
    a41 = 1 - (a42 + a43)
    a31 =
        (1 - 9/2 * γ + 3/2 * γ^2) * a42 +
        (11/4 - 21/2 * γ + 15/4 * γ^2) * a43 -
        7/2 + 13γ - 9/2 * γ^2
    a32 = cs[3] - a31
    A = SMatrix{4, 4}(
        # Column 1
        0, γ, a31, a41,
        # Column 2
        0, 0, a32, a42,
        # Column 3
        0, 0, 0, a43,
        # Column 4
        0, 0, 0, 0,
    )
    (; A, bs, cs,)
end

# Note: this is an ESDIRK scheme (explicit, singly diagonally implicit Runge-Kutta)
function tableau_implicit(::Ascher343)
    γ = 0.435866521508459
    b1 = -3/2 * γ^2 + 4γ - 1/4
    b2 = 1 - (γ + b1)
    cs = SVector(0.0, γ, (1 + γ) / 2, 1.0)
    bs = SVector(0, b1, b2, γ)
    A = SMatrix{4, 4}(
        # Column 1
        0, 0, 0, 0,
        # Column 2
        0, γ, (1 - γ) / 2, b1,
        # Column 3
        0, 0, γ, b2,
        # Column 4
        0, 0, 0, γ,
    )
    (; A, bs, cs,)
end
