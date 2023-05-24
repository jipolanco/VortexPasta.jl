"""
    ExplicitTemporalScheme

Abstract type defining an explicit temporal scheme.
"""
abstract type ExplicitTemporalScheme end

# By default, schemes allow changing the timestep.
can_change_dt(::ExplicitTemporalScheme) = true

"""
    TemporalSchemeCache{Scheme <: ExplicitTemporalScheme}

Contains buffers needed by a temporal scheme.
"""
struct TemporalSchemeCache{
        Scheme <: ExplicitTemporalScheme,
        Nf, Nv,
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfArray{<:Vec3},
    }
    scheme :: Scheme
    fc     :: NTuple{Nf, Filaments}
    vc     :: NTuple{Nv, Velocities}
end

scheme(c::TemporalSchemeCache) = c.scheme
can_change_dt(c::TemporalSchemeCache) = can_change_dt(scheme(c))

function init_cache(
        scheme::ExplicitTemporalScheme,
        fs::VectorOfFilaments, vs::VectorOfArray,
    )
    Nf = nbuf_filaments(scheme)
    Nv = nbuf_velocities(scheme)
    fc = ntuple(_ -> map(similar, fs), Val(Nf))
    vc = ntuple(_ -> similar(vs), Val(Nv))
    TemporalSchemeCache{
        typeof(scheme), Nf, Nv,
        typeof(fs), typeof(vs),  # needed in case tuples are empty (case of Euler)
    }(scheme, fc, vc)
end

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

function update_velocities!(
        rhs!::F, advect!::G, cache::TemporalSchemeCache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    resize!(cache, iter.fs)  # in case the number of nodes (or filaments) has changed
    _update_velocities!(scheme(cache), rhs!, advect!, cache, iter)
end

include("Euler.jl")
include("RK4.jl")
include("SSPRK33.jl")
include("DP5.jl")
