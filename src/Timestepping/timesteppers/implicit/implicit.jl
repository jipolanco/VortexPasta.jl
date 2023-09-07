"""
    ImplicitScheme <: TemporalScheme

Abstract type representing an implicit (possibly multi-stage) scheme.

These kinds of schemes are mainly meant to be used as inner schemes (for resolving
fast-evolving dynamics which are also cheap to compute) when using
[`MultirateSchemes`](@ref).
"""
abstract type ImplicitScheme <: TemporalScheme end

# Solve nonlinear system using fixed-point iterations (kind of brute-force method).
# We stop when the difference compared to the previous positions converges to a
# constant.
function solve_fixed_point!(
        ftmp, rhs!::F, advect!::G, iter::AbstractSolver, vtmp, v_explicit;
        dt = iter.time.dt, t = iter.time.t, cdt = dt, fbase, aI_diag,
        nmax = 200, rtol = 1e-10,
    ) where {F <: Function, G <: Function}
    @assert fbase !== ftmp
    vdiff_prev = vector_difference(fbase, ftmp)
    n = 0
    while n < nmax
        n += 1
        rhs!(vtmp, ftmp, t + cdt, iter)  # compute RHS at the latest location
        @. vtmp = v_explicit + aI_diag * vtmp  # full velocity estimate
        # Update guess for filament location
        advect!(ftmp, vtmp, dt; fbase)
        vdiff = vector_difference(fbase, ftmp)
        rdiff = abs(vdiff - vdiff_prev) / vdiff
        rdiff < rtol && break
        vdiff_prev = vdiff
    end
    ftmp
end

include("CrankNicolson.jl")
