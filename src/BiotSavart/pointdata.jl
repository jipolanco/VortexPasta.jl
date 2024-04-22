# This stores the locations s⃗ and charges q * s⃗′ used to compute both short-range and
# long-range interactions, allowing to reuse computations. Note that locations and charges
# are obtained via interpolation in-between filament nodes.
#
# This is also reused by long-range interactions to perform interpolations from Fourier to
# physical space (see `interpolate_to_physical!`).
# In that case, `points` contains the interpolation points (usually the filament nodes) and
# `charges` the interpolation values (usually velocities or streamfunction values).
struct PointData{
        T <: AbstractFloat,
        S <: Union{T, Complex{T}},
        Points <: StructVector{Vec3{T}},
        # Note: complex is needed by some long-range backends such as FINUFFT (even though values are always real!)
        Charges <: StructVector{Vec3{S}},
        Filament <: AbstractFilament,
        Segments <: AbstractVector{Segment{Filament}},
    }
    points   :: Points    # interpolated locations s⃗ on segments
    charges  :: Charges   # rescaled tangent vector q * s⃗′ on segments (where `q` is the quadrature weight)
    segments :: Segments  # filament segment on which each location s⃗ is located
end

function PointData(::Type{T}, ::Type{S}, ::Type{F}) where {T, S, F <: AbstractFilament}
    points = StructVector{Vec3{T}}(undef, 0)
    charges = StructVector{Vec3{S}}(undef, 0)
    Seg = Segment{F}
    @assert isconcretetype(Seg)
    segments = Seg[]
    PointData(points, charges, segments)
end

function Base.copy(data::PointData)
    (; points, charges, segments,) = data
    PointData(copy(points), copy(charges), copy(segments))
end

"""
    set_num_points!(data::PointData, N::Integer)

Set the total number of non-uniform points that the cache must hold.

This will reallocate space to make all points fit in the cache. It will also reset the
contributions of previously-added charges.

Must be called before [`add_point!`](@ref).
"""
function set_num_points!(data::PointData, N)
    resize!(data.points, N)
    resize!(data.charges, N)
    resize!(data.segments, N)
    data
end

function _count_charges(quad::StaticSizeQuadrature, fs::AbstractVector{<:ClosedFilament})
    Nq = length(quad)        # number of evaluation points per filament segment
    Np = sum(f -> length(segments(f)), fs; init = 0)  # total number of segments among all filaments (assumes closed filaments!!)
    Np * Nq
end

"""
    add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, quad::StaticSizeQuadrature)
    add_point_charges!(cache::LongRangeCache, fs::AbstractVector{<:AbstractFilament})

Add vector charges at multiple non-uniform locations.

This can be used for both short-range and long-range computations.

In the case of long-range computations, this must be done before type-1 NUFFTs, to transform
from non-uniform data in physical space to uniform data in Fourier space. It must be called
before [`compute_vorticity_fourier!`](@ref).
"""
function add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, quad::StaticSizeQuadrature)
    Ncharges = _count_charges(quad, fs)
    set_num_points!(data, Ncharges)
    inds = ChunkSplitters.chunks(fs; n = Threads.nthreads())
    Threads.@threads :static for subinds ∈ inds
        prev_indices = firstindex(fs):(first(subinds) - 1)  # filament indices given to all previous chunks
        n = _count_charges(quad, view(fs, prev_indices))    # we will start writing at index n + 1
        for i ∈ subinds
            n = _add_point_charges!(data, fs[i], n, quad)
        end
    end
    nothing
end

function _add_point_charges!(data::PointData, f, n::Int, quad::StaticSizeQuadrature)
    @assert eachindex(data.points) == eachindex(data.charges) == eachindex(data.segments)
    nlast = n + length(segments(f)) * length(quad)
    checkbounds(data.points, nlast)
    ζs, ws = quadrature(quad)
    ts = knots(f)
    @inbounds for (i, seg) ∈ pairs(segments(f))
        Δt = ts[i + 1] - ts[i]
        for (ζ, w) ∈ zip(ζs, ws)
            s⃗ = f(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
            # Note: the vortex circulation Γ is included in the Ewald operator and
            # doesn't need to be included here.
            q = w * Δt
            add_pointcharge!(data, s⃗, q * s⃗′, seg, n += 1)
        end
    end
    @assert n == nlast
    n
end

function _add_point_charges!(data::PointData, f, n::Int, ::NoQuadrature)
    @assert eachindex(data.points) == eachindex(data.charges) == eachindex(data.segments)
    nlast = n + length(segments(f))
    checkbounds(data.points, nlast)
    @inbounds for (i, seg) ∈ pairs(segments(f))
        s⃗ = (f[i] + f[i + 1]) ./ 2
        s⃗′_dt = f[i + 1] - f[i]
        add_pointcharge!(data, s⃗, s⃗′_dt, seg, n += 1)
    end
    @assert n == nlast
    n
end

function add_pointcharge!(data::PointData, X::Vec3, Q::Vec3, s::Segment, i::Int)
    @inbounds data.points[i] = X
    @inbounds data.charges[i] = Q
    @inbounds data.segments[i] = s
    data
end

"""
    add_point!(data::PointData, X::Vec3, i::Int)

Add an interpolation point for type-2 NUFFT.

This is used for type-2 NUFFTs, to transform (interpolate) from uniform data in Fourier
space to non-uniform data in physical space.

The total number of locations must be first set via [`set_num_points!`](@ref).
"""
function add_point!(data::PointData, X::Vec3, i::Int)
    @inbounds data.points[i] = X
    data
end
