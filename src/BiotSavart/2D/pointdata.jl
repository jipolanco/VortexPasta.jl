struct PointVortexData{
        T <: AbstractFloat,
        Vecs <: StructVector{Vec2{T}},
        Signs <: AbstractVector{SignValue},
        HostBuf <: HostVector{T},
    } <: AbstractPointData{T}
    points::Vecs       # vortex locations in the actual physical domain of periods Ls = (Lx, Ly)
    points_mod::Vecs   # vortex locations to be passed to the long-range backend (in [0, 2π]² for NonuniformFFTsBackend)
    charges::Signs
    buf_host::HostBuf  # buffer for host-device transfers
end

get_source_points_shortrange(data::PointVortexData) = data.points       # input points for short-range backend (not rescaled)
get_destination_points_shortrange(data::PointVortexData) = data.points  # output points for short-range backend (not rescaled)

get_source_points_longrange(data::PointVortexData) = data.points_mod      # input points for long-range backend (already rescaled if necessary)
get_destination_points_longrange(data::PointVortexData) = data.points_mod  # output points for long-range backend (already rescaled if necessary)

@inline function Adapt.adapt_structure(to::KA.Backend, p::PointVortexData)
    PointVortexData(
        adapt(to, p.points),
        adapt(to, p.points_mod),
        adapt(to, p.charges),
        adapt(to, p.buf_host),
    )
end

function PointVortexData(::Type{T}) where {T <: AbstractFloat}
    points = StructVector{Vec2{T}}(undef, 0)
    points_mod = similar(points)
    charges = Vector{SignValue}(undef, 0)
    buf_host = HostVector{T}(undef, CPU(), 0)  # not very useful to have this with backend = CPU, but adapt_structure will take care of changing the backend
    PointVortexData(points, points_mod, charges, buf_host)
end

function Base.copy(data::PointVortexData)
    buf_host = similar(data.buf_host, 0)
    PointVortexData(copy(data.points), copy(data.points_mod), copy(data.charges), buf_host)
end

# This may be used by the long-range backend to rescale coordinates to match the backend's
# expected domain period (2π in the case of NonuniformFFTs).
function _rescale_coordinates!(data::PointVortexData, Ls::NTuple, L_expected::Real)
    (; points, points_mod) = data
    for (xs, xs_out, L) in zip(StructArrays.components(points), StructArrays.components(points_mod), Ls)
        # We always keep a copy of the original nodes.
        resize_no_copy!(xs_out, length(xs))
        _rescale_coordinates!(xs_out, xs, L, L_expected)
    end
    nothing
end

# Backend doesn't define `expected_period`, so no rescaling is needed.
# We only copy points -> points_mod.
# This is typically the case of ExactSumBackend.
function _rescale_coordinates!(data::PointVortexData, Ls::NTuple, L_expected::Nothing)
    (; points, points_mod) = data
    resize_no_copy!(points_mod, length(points))
    copyto!(points_mod, points)
    nothing
end
