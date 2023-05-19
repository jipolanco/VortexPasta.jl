abstract type ExplicitTemporalScheme end

abstract type TemporalSchemeCache end

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
