abstract type TemporalScheme end

# By default, schemes allow changing the timestep.
can_change_dt(::TemporalScheme) = true

"""
    TemporalSchemeCache{Scheme <: TemporalScheme}

Contains buffers needed by a temporal scheme.
"""
struct TemporalSchemeCache{
        Scheme <: TemporalScheme,
        Nf, Nv,
        Filaments <: VectorOfFilaments,
        Velocities <: VectorOfVectors{<:Vec3},
    }
    scheme :: Scheme
    fc     :: NTuple{Nf, Filaments}
    vc     :: NTuple{Nv, Velocities}
end

scheme(c::TemporalSchemeCache) = c.scheme
can_change_dt(c::TemporalSchemeCache) = can_change_dt(scheme(c))

function init_cache(
        scheme::TemporalScheme,
        fs::VectorOfFilaments, vs::VectorOfVectors,
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

    # 1. Resize vectors of vectors if number of filaments changed.
    map(buf -> resize_container!(buf, fs), fc)
    map(buf -> resize_container!(buf, fs), vc)

    # 2. Resize individual vectors if number of nodes in each filament changed.
    # This is similar to resize_contained_vectors!.
    for (i, f) ∈ pairs(fs)
        N = length(f)
        for fbuf ∈ fc
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

include("explicit/explicit.jl")
include("imex/imex.jl")  # implicit-explicit
