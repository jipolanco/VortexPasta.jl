"""
    PointData

Stores point data (values on filaments) used to compute short-range and long-range interactions.

This includes the locations `s⃗` and charges `q * s⃗′` used to compute both short-range and
long-range interactions, allowing to reuse computations. Note that locations and charges
are obtained via interpolation in-between filament nodes.

This is also reused by long-range interactions to perform interpolations from Fourier to
physical space (see [`interpolate_to_physical!`](@ref)). In that case, `points` contains the
interpolation points (usually the filament nodes) and `charges` the interpolation values
(usually velocities or streamfunction values).

# Construction

    PointData(::Type{T}, ::Type{S}, ::Type{F}) -> PointData

where:

- `T` is a float type (usually `Float32` or `Float64`).

- `S` is either `T` (for real-valued charges) or `Complex{T}` (for complex-valued charges).
  Usually `T` is fine; `Complex{T}` may be needed by specific long-range backends (such as
  `FINUFFTBackend`, which is no longer available).

- `F <: AbstractFilament` is the filament type (e.g. `ClosedFilament{…}`).
"""
struct PointData{
        T <: AbstractFloat,
        S <: Union{T, Complex{T}},
        Points <: StructVector{Vec3{T}},
        # Note: complex is needed by some long-range backends such as FINUFFT (even though values are always real!)
        Charges <: StructVector{Vec3{S}},
        PointsHost <: StructVector{Vec3{T}},
        ChargesHost <: StructVector{Vec3{S}},
        Filament <: AbstractFilament,
        Segments <: AbstractVector{Segment{Filament}},
    }
    points    :: Points      # interpolated locations s⃗ on segments
    charges   :: Charges     # rescaled tangent vector q * s⃗′ on segments (where `q` is the quadrature weight)
    points_h  :: PointsHost  # CPU buffer which may be used for intermediate host-device transfers
    charges_h :: ChargesHost # CPU buffer which may be used for intermediate host-device transfers
    segments  :: Segments    # filament segment on which each location s⃗ is located
end

# If `to` corresponds to a GPU backend, create PointData object on the GPU (useful for long-range computations
# with GPU backend). Returns `p` without allocations if data must be adapted from CPU to
# CPU. See https://cuda.juliagpu.org/dev/tutorials/custom_structs/.
@inline function Adapt.adapt_structure(to, p::PointData)
    PointData(
        adapt(to, p.points),
        adapt(to, p.charges),
        p.points_h,   # this is always on the CPU
        p.charges_h,  # this is always on the CPU
        p.segments,   # for now, keep segments on the CPU (we don't use them on the GPU)
    )
end

function PointData(::Type{T}, ::Type{S}, ::Type{F}) where {T, S, F <: AbstractFilament}
    points = StructVector{Vec3{T}}(undef, 0)
    charges = StructVector{Vec3{S}}(undef, 0)
    Seg = Segment{F}
    @assert isconcretetype(Seg)
    segments = Seg[]
    PointData(points, charges, copy(points), copy(charges), segments)
end

function Base.copy(data::PointData)
    (; points, charges, segments,) = data
    points_h = similar(data.points_h, 0)   # empty arrays (we don't need them to be identical to the original ones)
    charges_h = similar(data.charges_h, 0)
    PointData(copy(points), copy(charges), points_h, charges_h, copy(segments))
end

# This is useful in particular for host -> device copies.
# Note that arrays are resized to match those in `src`.
function Base.copy!(dst::PointData, src::PointData)
    copy!(dst.points, src.points)
    copy!(dst.charges, src.charges)
    # Note that both `segments` fields may point to the same object; see `adapt_structure` above.
    dst.segments === src.segments || copy!(dst.segments, src.segments)  # not useful on the GPU...
    dst
end

"""
    set_num_points!(data::PointData, N::Integer)

Set the total number of non-uniform points that the cache must hold.

This will reallocate space to make all points fit in the cache. It will also reset the
contributions of previously-added charges.
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
    inds = OhMyThreads.chunks(eachindex(fs); n = Threads.nthreads())
    Threads.@threads for subinds ∈ inds
        subinds === nothing && continue  # subinds can be nothing; see iterate(c::Chunk) in ChunkSplitters
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
