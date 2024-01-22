# Write data on filament nodes -- case of scalar data (e.g. local curvature magnitude)
function write_data_on_filaments(
        ::Type{T}, gdata, vs::AbstractVector{<:AbstractVector}, name, refinement,
    ) where {T <: Number}
    @assert T === eltype(eltype(vs))
    num_points = sum(v -> refinement * length(v) + 1, vs)  # the +1 is to include the endpoint
    dtype = HDF5.datatype(T)
    dspace = HDF5.dataspace((num_points,))
    dset = HDF5.create_dataset(gdata, name, dtype, dspace)
    info = (; dtype, dspace, dset,)
    ldims = ()  # scalar data: 0-dimensional elements

    n = if refinement == 1
        write_data_unrefined(info, vs; ldims)
    else
        write_data_refined(info, vs, refinement; ldims)
    end
    @assert n == num_points

    map(close, info)

    nothing
end

# Write data on filament nodes -- case of vector data (e.g. velocity)
function write_data_on_filaments(
        ::Type{V}, gdata, vs::AbstractVector{<:AbstractVector}, name, refinement,
    ) where {V <: AbstractVector}
    @assert V === eltype(eltype(vs))
    T = eltype(V)
    @assert T <: Number
    N = length(V)  # usually == 3 (works when V <: StaticArray)
    num_points = sum(v -> refinement * length(v) + 1, vs)  # the +1 is to include the endpoint
    dtype = HDF5.datatype(T)
    dspace = HDF5.dataspace((N, num_points))
    dset = HDF5.create_dataset(gdata, name, dtype, dspace)
    info = (; dtype, dspace, dset,)
    ldims = (N,)   # vector data

    n = if refinement == 1
        write_data_unrefined(info, vs; ldims)
    else
        write_data_refined(info, vs, refinement; ldims)
    end

    @assert n == num_points
    map(close, info)

    nothing
end

function write_data_unrefined(info, vs; ldims::Dims = ())
    (; dtype, dspace, dset,) = info
    memtype = dtype
    memspace_node = HDF5.dataspace(ldims)  # this is for writing a value on a single point
    lranges = map(N -> 1:N, ldims)
    n = 0
    for v ∈ vs
        Np = length(v)
        # Write all data associated to a single filament at once
        memspace_filament = HDF5.dataspace((ldims..., Np,))
        HDF5.select_hyperslab!(dspace, (lranges..., (n + 1):(n + Np),))
        HDF5.API.h5d_write(dset.id, memtype.id, memspace_filament.id, dspace.id, dset.xfer, v)
        close(memspace_filament)
        # Close the loop: this writes v[begin]
        n += Np + 1
        HDF5.select_hyperslab!(dspace, (lranges..., n,))
        HDF5.API.h5d_write(dset.id, memtype.id, memspace_node.id, dspace.id, dset.xfer, v)
    end
    close(memspace_node)
    n
end

function write_data_refined(info, vs, refinement; ldims::Dims = ())
    V = eltype(eltype(vs))  # usually SVector{3, T} or just T (scalar data)
    T = eltype(V)
    @assert T <: Number
    n = 0
    buf = Bumper.default_buffer()
    for v ∈ vs
        @no_escape buf begin
            Np = refinement * length(v) + 1  # number of output points (the +1 is to close the loop)
            Xs = @alloc(T, ldims..., Np)
            refine_data_on_filament!(Xs, v, refinement)
            write_points_to_hdf5(info, Xs, n, ldims, Np)
            n += Np
        end
    end
    n
end

function refine_data_on_filament!(Xs::AbstractArray, v::AbstractVector, refinement::Int)
    n = firstindex(Xs)::Int - 1
    for j ∈ eachindex(v), m ∈ 1:refinement
        u = _interpolate_on_segment(v, j, m, refinement)
        for w ∈ u
            Xs[n += 1] = w
        end
    end
    for w ∈ v[end]
        Xs[n += 1] = w  # close the loop
    end
    @assert n == lastindex(Xs)
    Xs
end

# Interpolate value on segment. For now we simply use linear interpolation (it's just for
# visualisation anyways).
function _interpolate_on_segment(v, i, m, M)
    α = (m - 1) / M
    (1 - α) * v[i] + α * _vnext(v, i)
end

_vnext(v::PaddedVector, i) = v[i + 1]
_vnext(v::AbstractVector, i) = v[i == lastindex(v) ? firstindex(v) : (i + 1)]
