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
        fs::AbstractVector{<:AbstractFilament},
        gs::AbstractVector{<:AbstractFilament},
    )
    sqdiff = 0.0
    for (f, g) ∈ zip(fs, gs)
        for i ∈ eachindex(f, g)
            @inbounds sqdiff += sum(abs2, f[i] - g[i])
        end
    end
    sqrt(sqdiff)
end

include("Euler.jl")
