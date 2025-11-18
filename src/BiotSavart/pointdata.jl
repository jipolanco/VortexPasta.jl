"""
    PointData

Stores point data (values on filaments) used to compute short-range and long-range interactions.

Among the stored fields are:

- `nodes::StructVector{Vec3{T}}`: filament nodes (discretisation points) where fields like
  velocity and streamfunction is to be computed;

- `points::StructVector{Vec3{T}}`: quadrature points on filament segments where "vorticity"
  is to be evaluated;

- `charges::StructVector{Vec3{T}}`: "vorticity" vector evaluated on quadrature points.
  More precisely, this is `w * s⃗′` where `w` is a quadrature weight and `s⃗′` is the local
  tangent vector.

Quadrature points and charges are obtained via interpolation in-between filament nodes.

# Construction

    PointData(::Type{T}, ::Type{F}) -> PointData

where:

- `T <: AbstractFloat` is a float type (usually `Float32` or `Float64`).

- `F <: AbstractFilament` is the filament type (e.g. `ClosedFilament{…}`).
"""
struct PointData{
        T <: AbstractFloat,
        Vecs <: StructVector{Vec3{T}},
        VecsHost <: StructVector{Vec3{T}},
        Filament <: AbstractFilament,
        Segments <: AbstractVector{Segment{Filament}},
    }
    nodes     :: Vecs      # filament nodes (where velocity will be computed)
    points    :: Vecs      # interpolated locations s⃗ on segments
    charges   :: Vecs      # rescaled tangent vector q * s⃗′ on segments (where `q` is the quadrature weight)
    points_h  :: VecsHost  # CPU buffer which may be used for intermediate host-device transfers
    charges_h :: VecsHost  # CPU buffer which may be used for intermediate host-device transfers
    segments  :: Segments  # filament segment on which each location s⃗ is located
end

# If `to` corresponds to a GPU backend, create PointData object on the GPU (useful for long-range computations
# with GPU backend). Returns `p` without allocations if data must be adapted from CPU to
# CPU. See https://cuda.juliagpu.org/dev/tutorials/custom_structs/.
@inline function Adapt.adapt_structure(to, p::PointData)
    PointData(
        adapt(to, p.nodes),
        adapt(to, p.points),
        adapt(to, p.charges),
        p.points_h,   # this is always on the CPU
        p.charges_h,  # this is always on the CPU
        p.segments,   # for now, keep segments on the CPU (we don't use them on the GPU)
    )
end

function PointData(::Type{T}, ::Type{F}) where {T <: AbstractFloat, F <: AbstractFilament}
    nodes = StructVector{Vec3{T}}(undef, 0)
    points = similar(nodes)
    charges = similar(nodes)
    Seg = Segment{F}
    @assert isconcretetype(Seg)
    segments = Seg[]
    PointData(nodes, points, charges, copy(points), copy(charges), segments)
end

function Base.copy(data::PointData)
    (; nodes, points, charges, segments,) = data
    points_h = similar(data.points_h, 0)   # empty arrays (we don't need them to be identical to the original ones)
    charges_h = similar(data.charges_h, 0)
    PointData(copy(nodes), copy(points), copy(charges), points_h, charges_h, copy(segments))
end

# This is useful in particular for host -> device copies.
# Note that arrays are resized to match those in `src`.
function Base.copy!(dst::PointData, src::PointData)
    copy!(dst.nodes, src.nodes)
    copy!(dst.points, src.points)
    copy!(dst.charges, src.charges)
    # Note that both `segments` fields may point to the same object; see `adapt_structure` above.
    dst.segments === src.segments || copy!(dst.segments, src.segments)  # not useful on the GPU...
    dst
end

"""
    set_num_points!(data::PointData, Np::Integer, quad::StaticSizeQuadrature)

Set the total number of non-uniform points that the cache must hold.

This will reallocate space to make all points fit in the cache. It will also reset the
contributions of previously-added charges.
"""
function set_num_points!(data::PointData, Np, quad::StaticSizeQuadrature)
    Nq = Np * length(quad)  # number of quadrature nodes
    resize!(data.nodes, Np)
    resize!(data.points, Nq)
    resize!(data.charges, Nq)
    resize!(data.segments, Nq)
    data
end

# Total number of independent nodes among all filaments
_count_nodes(fs::AbstractVector{<:AbstractFilament}) = sum(f -> length(nodes(f)), fs; init = 0)

"""
    add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, Ls::NTuple, quad::StaticSizeQuadrature)
    add_point_charges!(cache::LongRangeCache, fs::AbstractVector{<:AbstractFilament})

Add vector charges at multiple non-uniform locations.

This can be used for both short-range and long-range computations.

In periodic domains, each point is folded into the main periodic cell (in ``[0, L]^3``).
This can be accounted for to speed-up operations done later with point data, in particular in
short-range computations when determining the minimal distance between two points in the
periodic lattice.

In the case of long-range computations, this must be done before type-1 NUFFTs, to transform
from non-uniform data in physical space to uniform data in Fourier space. It must be called
before [`compute_vorticity_fourier!`](@ref).
"""
function add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, Ls::NTuple, quad::StaticSizeQuadrature)
    Np = _count_nodes(fs)
    set_num_points!(data, Np, quad)
    chunks = FilamentChunkIterator(fs)  # defined in shortrange/shortrange.jl
    @sync for chunk in chunks
        isempty(chunk) && continue  # don't spawn a task if it will do no work
        Threads.@spawn let
            prev_indices = firstindex(fs):(first(chunk) - 1)  # filament indices given to all previous chunks
            n = _count_nodes(view(fs, prev_indices))  # we will start writing at index n + 1
            for i in chunk
                n = _add_point_charges!(data, fs[i], Ls, n, quad)
            end
        end
    end
    nothing
end

function _add_point_charges!(data::PointData, f, Ls, n::Int, quad::StaticSizeQuadrature)
    @assert eachindex(data.points) == eachindex(data.charges) == eachindex(data.segments)
    m = length(quad) * n  # current index in points/charges/segments vectors
    nlast = n + length(segments(f))
    mlast = length(quad) * nlast
    checkbounds(data.nodes, nlast)
    checkbounds(data.points, mlast)
    ζs, ws = quadrature(quad)
    ts = knots(f)
    @inbounds for (i, seg) ∈ pairs(segments(f))
        data.nodes[n += 1] = Filaments.fold_coordinates_periodic(f[i], Ls)
        Δt = ts[i + 1] - ts[i]
        for (ζ, w) ∈ zip(ζs, ws)
            s⃗ = f(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
            # Note: the vortex circulation Γ is included in the Ewald operator and
            # doesn't need to be included here.
            q = w * Δt
            _add_pointcharge!(data, s⃗, q * s⃗′, seg, Ls, m += 1)
        end
    end
    @assert n == nlast
    @assert m == mlast
    n
end

function _add_point_charges!(data::PointData, f, Ls, n::Int, ::NoQuadrature)
    @assert eachindex(data.nodes) == eachindex(data.points) == eachindex(data.charges) == eachindex(data.segments)
    nlast = n + length(segments(f))
    checkbounds(data.nodes, nlast)
    checkbounds(data.points, nlast)
    @inbounds for (i, seg) ∈ pairs(segments(f))
        n += 1
        data.nodes[n] = Filaments.fold_coordinates_periodic(f[i], Ls)
        s⃗ = (f[i] + f[i + 1]) ./ 2   # segment midpoint as quadrature node
        s⃗′_dt = f[i + 1] - f[i]
        _add_pointcharge!(data, s⃗, s⃗′_dt, seg, Ls, n)
    end
    @assert n == nlast
    n
end

function _add_pointcharge!(data::PointData, X::Vec3, Q::Vec3, s::Segment, Ls::NTuple{3}, i::Int)
    @inbounds data.points[i] = Filaments.fold_coordinates_periodic(X, Ls)
    @inbounds data.charges[i] = Q
    @inbounds data.segments[i] = s
    data
end
