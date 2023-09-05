"""
    MultirateScheme <: TemporalScheme

Abstract type defining a multirate scheme.

The idea is to treat different terms of the evolution equations with different timesteps
(and possibly different integration schemes). Concretely, the fast dynamics is integrated
with a smaller timestep than the slowly-evolving motions.
"""
abstract type MultirateScheme <: TemporalScheme end

# Adds "forcing" term (from slow dynamics) to the RHS of the ODE for the fast dynamics.
# The idea is to sum all terms by looping over all values just once. We do this using
# broadcasting, by first building a lazy Broadcasted object which can be then evaluated.
function _MRI_inner_add_forcing_term!(
        vs, vslow::NTuple, Γs::NTuple, τ,
    )
    i = length(vslow)  # this usually corresponds to the current outer RK stage
    bc_slow = ntuple(Val(i)) do j
        @inline
        coefs = map(Γ -> Γ[i, j], Γs)
        c = evalpoly(τ, coefs)  # = Γs[1][i, j] + Γs[2][i, j] * τ + ...
        Base.broadcasted(*, c, vslow[j])
    end
    bc = Base.broadcasted(+, vs, bc_slow...)  # sum all terms including the fast velocity `vs`
    Base.materialize!(vs, bc)  # evaluate the whole thing
    vs
end

include("midpoint.jl")
include("MRI-GARK-ERK33a.jl")
include("MRI-GARK-ERK45a.jl")
