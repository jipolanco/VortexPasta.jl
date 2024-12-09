using StaticArrays: SVector, SMatrix  # may be used for RK tableaus

"""
    TemporalScheme

Abstract type representing a timestepping scheme.
"""
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
        Filaments <: VectorOfVectors{<:Vec3, <:AbstractFilament},
        Velocities <: VectorOfVectors{<:Vec3},
    }
    scheme :: Scheme
    fc     :: NTuple{Nf, Filaments}
    vc     :: NTuple{Nv, Velocities}

    function TemporalSchemeCache(scheme, fc::NTuple{Nf}, vc::NTuple{Nv}) where {Nf, Nv}
        @assert nbuf_filaments(scheme) == Nf
        @assert nbuf_velocities(scheme) == Nv
        new{
            typeof(scheme), Nf, Nv, eltype(fc), eltype(vc),
        }(scheme, fc, vc)
    end
end

Base.summary(io::IO, c::TemporalSchemeCache) = print(io, "TemporalSchemeCache(", scheme(c), ")")

scheme(c::TemporalSchemeCache) = c.scheme
can_change_dt(c::TemporalSchemeCache) = can_change_dt(scheme(c))

function init_cache(
        scheme::TemporalScheme,
        fs::VectorOfVectors, vs::VectorOfVectors,
    )
    Nf = nbuf_filaments(scheme)
    Nv = nbuf_velocities(scheme)
    fc = ntuple(_ -> similar(fs), Val(Nf))
    vc = ntuple(_ -> similar(vs), Val(Nv))
    TemporalSchemeCache(scheme, fc, vc)
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
        vL::VectorOfVectors, rhs!::F, advect!::G, cache::TemporalSchemeCache, iter::AbstractSolver;
        resize_cache = true, kws...,
    ) where {F <: Function, G <: Function}
    if resize_cache
        resize!(cache, iter.fs)  # in case the number of nodes (or filaments) has changed
    end
    _update_velocities!(scheme(cache), vL, rhs!, advect!, cache, iter; kws...)
end

include("explicit/explicit.jl")
include("implicit/implicit.jl")
include("imex/imex.jl")  # implicit-explicit
include("multirate/multirate.jl")
include("split/split.jl")  # e.g. Strang splitting
