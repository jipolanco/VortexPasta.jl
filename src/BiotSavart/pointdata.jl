"""
    PointData{T <: AbstractFloat}
    PointData(::Type{T}) -> PointData{T}

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
"""
struct PointData{
        T <: AbstractFloat,
        Vecs <: StructVector{Vec3{T}},
        Scalars <: AbstractVector{T},
        Indices <: AbstractVector{<:Integer},
        VecsHost <: StructVector{Vec3{T}},
    }
    nodes         :: Vecs     # [Np] filament nodes (where velocity will be computed)
    node_idx_prev :: Indices  # [Np] index in 1:Np of node located to the "left" of a given node (accounts for periodic wrapping of closed filaments) -> only used in short-range when avoid_explicit_erf == false
    points    :: Vecs      # [Nq] quadrature points, i.e. interpolated locations s⃗ on segments
    charges   :: Vecs      # [Nq] rescaled tangent vector q * s⃗′ on segments (where `q` is the quadrature weight)
    derivatives_on_nodes :: NTuple{2, Vecs}   # [Np] derivatives on filament nodes (s′, s″)
    subsegment_lengths :: NTuple{2, Scalars}  # [Np] lengths δ⁻ and δ⁺ of subsegments for local BS term (accounting for lia_segment_fraction)
    # TODO: do we need this?:
    charges_h :: VecsHost  # CPU buffer which may be used for intermediate host-device transfers
end

# If `to` corresponds to a GPU backend, create PointData object on the GPU (useful for long-range computations
# with GPU backend). Returns `p` without allocations if data must be adapted from CPU to
# CPU. See https://cuda.juliagpu.org/dev/tutorials/custom_structs/.
@inline function Adapt.adapt_structure(to, p::PointData)
    PointData(
        adapt(to, p.nodes),
        adapt(to, p.node_idx_prev),
        adapt(to, p.points),
        adapt(to, p.charges),
        adapt(to, p.derivatives_on_nodes),
        adapt(to, p.subsegment_lengths),
        p.charges_h,  # this is always on the CPU
    )
end

function PointData(::Type{T}) where {T <: AbstractFloat}
    nodes = StructVector{Vec3{T}}(undef, 0)
    node_idx_prev = similar(nodes, Int32)  # Int32 should be enough for all practical purposes
    points = similar(nodes)
    charges = similar(nodes)
    derivatives_on_nodes = ntuple(_ -> similar(nodes), Val(2))
    subsegment_lengths = ntuple(_ -> similar(nodes, T), Val(2))
    PointData(nodes, node_idx_prev, points, charges, derivatives_on_nodes, subsegment_lengths, copy(charges))
end

function Base.copy(data::PointData)
    (; nodes, node_idx_prev, points, derivatives_on_nodes, subsegment_lengths, charges,) = data
    charges_h = similar(data.charges_h, 0)
    PointData(
        copy(nodes), copy(node_idx_prev), copy(points), copy(charges),
        map(copy, derivatives_on_nodes), map(copy, subsegment_lengths), charges_h
    )
end

