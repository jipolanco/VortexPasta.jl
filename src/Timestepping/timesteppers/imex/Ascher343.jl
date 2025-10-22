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
        scheme::Ascher343, vs, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs,) = iter
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
        @. v_explicit = (  # known part of the velocity
            AE[i, 1] * vE[1]
          # + AI[i, 1] * vI[1]  # note: AI[:, 1] == 0
        )
        @. vtmp = v_explicit + vI[i - 1] * AI[i, i]  # initial guess for effective velocity (use implicit value from previous stage)
        advect!(ftmp, vtmp, dt; fbase)  # initial guess for locations at stage i
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
        @. v_explicit = (  # known part of the velocity
            AE[i, 1] * vE[1]
          # + AI[i, 1] * vI[1]
          + AE[i, 2] * vE[2]
          + AI[i, 2] * vI[2]
        )
        @. vtmp = v_explicit + vI[i - 1] * AI[i, i]  # initial guess for effective velocity (use implicit value from previous stage)
        advect!(ftmp, vtmp, dt; fbase)  # initial guess for locations at stage i
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
        @. v_explicit = (  # known part of the velocity
            AE[i, 1] * vE[1]
          # + AI[i, 1] * vI[1]
          + AE[i, 2] * vE[2]
          + AI[i, 2] * vI[2]
          + AE[i, 3] * vE[3]
          + AI[i, 3] * vI[3]
        )
        @. vtmp = v_explicit + vI[i - 1] * AI[i, i]  # initial guess for effective velocity (use implicit value from previous stage)
        advect!(ftmp, vtmp, dt; fbase)  # initial guess for locations at stage i
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
