struct PointVortexData{
        T <: AbstractFloat,
        Vecs <: StructVector{Vec2{T}},
        Signs <: AbstractVector{SignValue},
        HostBuf <: HostVector{T},
       } <: AbstractPointData{T}
    points::Vecs
    charges::Signs
    buf_host::HostBuf  # buffer for host-device transfers
end

@inline function Adapt.adapt_structure(to::KA.Backend, p::PointVortexData)
    PointVortexData(
        adapt(to, p.points),
        adapt(to, p.charges),
        adapt(to, p.buf_host),
    )
end

function PointVortexData(::Type{T}) where {T <: AbstractFloat}
    points = StructVector{Vec2{T}}(undef, 0)
    charges = Vector{SignValue}(undef, 0)
    buf_host = HostVector{T}(undef, CPU(), 0)  # not very useful to have this with backend = CPU, but adapt_structure will take care of changing the backend
    PointVortexData(points, charges, buf_host)
end

function Base.copy(data::PointVortexData)
    buf_host = similar(data.buf_host, 0)
    PointVortexData(copy(data.points), copy(data.charges), buf_host)
end
