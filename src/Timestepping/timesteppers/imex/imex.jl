"""
    ImplicitExplicitScheme

Abstract type defining an implicit-explicit (a.k.a. IMEX) temporal scheme.

The defined IMEX schemes treat the localised induction approximation (LIA) term as implicit.
This may allow to increase the timestep, as it is the LIA term which imposes a small timestep
in VFM simulations. Moreover, since this term is quite cheap to compute (compared to
non-local interactions), treating it implicitly is not very expensive.
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
    @assert fbase !== ftmp
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

# Compute explicit part of the velocity at stage `i` of IMEX-RK scheme.
# This basically computes
#
#   v_explicit = sum(AE[i, j] * vE[j] + AI[i, j] * vI[j]) for j ∈ 1:(i - 1)
#
# in an efficient way, by constructing a lazy broadcast expression and then iterating over
# all arrays just once. In fact, it is the same as doing:
#
#   @. v_explicit = AE[i, 1] * vE[1] + AI[i, 1] * vI[i, 1] +
#                   AE[i, 2] * vE[2] + AI[i, 2] * vI[i, 2] +
#                   [...]                                  +
#                   AE[i, i - 1] * vE[i - 1] + AI[i, i - 1] * vI[i, i - 1]
#
function _imex_explicit_rhs!(::Val{i}, v_explicit, (AE, AI), (vE, vI)) where {i}
    bcs = ntuple(Val(i - 1)) do j
        Base.broadcasted(
            +,
            Base.broadcasted(*, AE[i, j], vE[j]),
            Base.broadcasted(*, AI[i, j], vI[j]),
        )
    end
    bc = Base.broadcasted(+, bcs...)
    Base.materialize!(v_explicit, bc)
    v_explicit
end

include("Euler.jl")
include("KenCarp3.jl")
include("KenCarp4.jl")
include("Ascher343.jl")
