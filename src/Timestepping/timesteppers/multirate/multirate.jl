"""
    MultirateScheme <: TemporalScheme

Abstract type defining a multirate scheme.

The idea is to treat different terms of the evolution equations with different timesteps
(and possibly different integration schemes). Concretely, the fast dynamics is integrated
with a smaller timestep (and a different *inner* scheme) than the slowly-evolving motions.
"""
abstract type MultirateScheme <: TemporalScheme end

inner_scheme(scheme::MultirateScheme) = scheme.inner

function Base.show(io::IO, scheme::MultirateScheme)
    print(io, nameof(typeof(scheme)), '(', scheme.inner, ", ", scheme.nsubsteps, ')')
end

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

# Run a single MRI-GARK stage.
# This corresponds to solving an auxiliary ODE in the range (ts[i], ts[i + 1]).
function _MRI_stage!(
        ::Val{i}, rhs!::F, advect!::G, ftmp, vS::NTuple, ts, iter,
        vs, cache_inner, Mfast, hfast, Γs,
    ) where {i, F, G}
    tsub = ts[i]
    cdt = ts[i + 1] - ts[i]
    vslow = ntuple(j -> vS[j], Val(i))
    function rhs_inner!(vs, fs, t, iter)
        rhs!(vs, fs, t, iter; component = Val(:fast))
        if i === 1
            @. vs = vs + vS[1]  # the first stage is always very simple
        else
            τ = (t - ts[i]) / cdt  # normalised time in [0, 1]
            _MRI_inner_add_forcing_term!(vs, vslow, Γs, τ)
        end
    end
    for _ ∈ 1:Mfast
        rhs_inner!(vs, ftmp, tsub, iter)
        update_velocities!(
            vs,
            rhs_inner!, advect!, cache_inner, iter;
            resize_cache = false, t = tsub, dt = hfast, fs = ftmp,
        )
        advect!(ftmp, vs, hfast; fbase = ftmp)
        tsub += hfast
    end
    @assert tsub ≈ ts[i + 1]
    nothing
end

include("midpoint.jl")
include("MRI-GARK-ERK33a.jl")
include("MRI-GARK-ERK45a.jl")
