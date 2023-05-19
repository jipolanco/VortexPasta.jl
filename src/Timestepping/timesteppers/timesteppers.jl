"""
    ExplicitTemporalScheme

Abstract type defining an explicit temporal scheme.
"""
abstract type ExplicitTemporalScheme end

# By default, schemes allow changing the timestep.
can_change_dt(::ExplicitTemporalScheme) = true

"""
    TemporalSchemeCache

Abstract type defining the cache of a temporal scheme.
"""
abstract type TemporalSchemeCache end

can_change_dt(c::TemporalSchemeCache) = can_change_dt(scheme(c))

function Base.resize!(cache::TemporalSchemeCache, fs::VectorOfFilaments)
    # TODO add support for changing the number of filaments...
    (; fc, vc,) = cache
    @assert length(fc) == length(fs)
    for (i, f) ∈ pairs(fs)
        N = length(f)
        resize!(fc[i], N)
        for v ∈ vc
            resize!(v[i], N)
        end
    end
    cache
end

include("RK4.jl")
