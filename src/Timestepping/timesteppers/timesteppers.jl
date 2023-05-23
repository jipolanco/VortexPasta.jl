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
    (; fc, vc,) = cache
    for (i, f) ∈ pairs(fs)
        N = length(f)
        for fbuf ∈ fc
            # TODO add support for changing the number of filaments...
            @assert length(fbuf) == length(fs)
            resize!(fbuf[i], N)
        end
        for vbuf ∈ vc
            resize!(vbuf[i], N)
        end
    end
    cache
end

include("Euler.jl")
include("RK4.jl")
