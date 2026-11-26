function read_point_data!(g::HDF5.Group, name, xs::AbstractVector)
    V = eltype(xs)  # usually SVector{3, T} or T <: Real
    _read_point_data!(V, g, name, xs)
end

# Scalar data or 3D vector data (easy)
# We use the low-level HDF5 API to avoid allocations.
_read_point_data!(::Type{T}, g::HDF5.Group, name, xs) where {T <: Number} =
    _read_full_dataset!(T, g, name, xs)

_read_point_data!(::Type{SVector{3, T}}, g::HDF5.Group, name, xs) where {T} =
    _read_full_dataset!(T, g, name, xs)

function _read_full_dataset!(::Type{T}, g, name, xs) where {T <: Number}
    V = eltype(xs)
    dims = if V <: SVector
        (length(V), size(xs)...)  # dimensions of HDF5 dataset
    else
        size(xs)
    end
    dset = HDF5.open_dataset(g, name)
    @assert size(dset) == dims
    memtype = HDF5.datatype(T)
    memspace = HDF5.dataspace(dims)
    filespace = memspace
    HDF5.API.h5d_read(dset.id, memtype.id, memspace.id, filespace.id, dset.xfer, xs)
    close(memtype)
    close(memspace)
    close(dset)
    xs
end

# 2D vector data: we need to discard the 3rd dimension
function _read_point_data!(::Type{SVector{2, T}}, g::HDF5.Group, name, xs) where {T}
    dims_mem = (2, size(xs)...)
    dims_file = (3, size(xs)...)  # dataset dimension in file
    dset = HDF5.open_dataset(g, name)
    @assert size(dset) == dims_file
    memtype = HDF5.datatype(T)
    memspace = HDF5.dataspace(dims_mem)
    filespace = HDF5.dataspace(dims_file)

    # Read hyperslab discarding 3rd dimension
    HDF5.select_hyperslab!(filespace, (1:2, eachindex(xs)))
    HDF5.API.h5d_read(dset.id, memtype.id, memspace.id, filespace.id, dset.xfer, xs)

    close(filespace)
    close(memspace)
    close(memtype)
    close(dset)
    xs
end
