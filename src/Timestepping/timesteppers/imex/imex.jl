"""
    ImplicitExplicitScheme <: TemporalScheme

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

function imex_rhs_implicit(rhs!::F) where {F}
    (args...) -> rhs!(args...; component = Val(:fast))
end

function imex_rhs_explicit(rhs!::F) where {F}
    (args...) -> rhs!(args...; component = Val(:slow))
end

function imex_rhs_full(rhs!::F) where {F}
    (args...) -> rhs!(args...; component = Val(:full))
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
