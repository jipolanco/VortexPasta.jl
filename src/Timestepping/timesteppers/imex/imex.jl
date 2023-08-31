"""
    ImplicitExplicitScheme

Abstract type defining an implicit-explicit (a.k.a. IMEX) temporal scheme.

The defined IMEX schemes treat the localised induction approximation (LIA) term as implicit.
This allows to increase the timestep, as it is the LIA term which imposes a small timestep
in fully explicit schemes. Moreover, since this term is quite cheap to compute (compared to
non-local interactions), treating it implicitly is not very expensive.

Note that, for now, IMEX schemes are only guaranteed to work when vortex filaments are
discretised using cubic splines.
"""
abstract type ImplicitExplicitScheme <: TemporalScheme end

function vector_difference(
        fs::AbstractVector{<:AbstractVector{<:SVector}},
        gs::AbstractVector{<:AbstractVector{<:SVector}},
    )
    sqdiff = 0.0
    for (f, g) ∈ zip(fs, gs)
        for i ∈ eachindex(f, g)
            @inbounds sqdiff += sum(abs2, f[i] - g[i])
        end
    end
    sqrt(sqdiff)
end

# Solve nonlinear system using fixed-point iterations (kind of brute-force method).
# We stop when the difference compared to the previous positions converges to a
# constant.
function solve_fixed_point!(
        ftmp, rhs!::F, advect!::G, iter::AbstractSolver, vtmp, v_explicit;
        cdt, fbase, aI_diag, nmax = 100, rtol = 1e-12,
    ) where {F <: Function, G <: Function}
    (; t, dt,) = iter.time
    vdiff_prev = vector_difference(fbase, ftmp)
    n = 0
    while n < nmax
        n += 1
        # Compute fast component at the latest location
        rhs!(vtmp, ftmp, t + cdt, iter; component = Val(:fast))
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

include("Euler.jl")
include("KenCarp3.jl")
include("KenCarp4.jl")
