export init_vtkhdf, write_point_data

using HDF5: HDF5

# This is to make sure we write the dataset type ("UnstructuredGrid") as ASCII.
# This is required by VTK/ParaView.
function datatype_ascii(s::AbstractString)
    dtype = HDF5.datatype(s)
    HDF5.API.h5t_set_cset(dtype.id, HDF5.API.H5T_CSET_ASCII)
    dtype
end

function open_or_create_group(parent, name)
    if haskey(parent, name)
        HDF5.open_group(parent, name)
    else
        HDF5.create_group(parent, name)
    end
end

"""
    init_vtkhdf(io::HDF5.File, fs::AbstractVector{<:AbstractFilament})

Initialise a new VTK HDF file with a list of filaments.

A VTK HDF file is an HDF5 file organised in such a way that it can be readily visualised in
tools such as ParaView.
Data in the file can be readily accessed using HDF5 tools and libraries.

Following the VTK HDF specification, this function creates a "VTKHDF" group on top of the
HDF5 file. Then, it creates the datasets allowing to describe the filaments as an
unstructured grid.

Some relevant datasets which are written are:

- `/VTKHDF/Points`: contains the coordinates of all filament nodes. Points are represented as
  an array of dimensions `(3, Np)` where `Np` is the total number of nodes. Note that the
  points *include the endpoint*, meaning that for a closed filament the initial coordinate
  appears twice in the file. This is done to disambiguate between *closed* and *infinite but
  unclosed* filaments (see [`end_to_end_offset`](@ref)).

- `/VTKHDF/Offsets`: contains the offset associated to each filament within the `Points`
  dataset. Its length is `1 + Nf` where `Nf` is the number of filaments. The points
  associated to the `i`-th filament are those with indices `(Offsets[i] + 1):Offsets[i + 1]`.
  Once again, note that number of nodes of filament `i` *including the endpoint* is
  `Offsets[i + 1] - Offsets[i]`.

After calling this function, one may use [`write_point_data`](@ref) to attach data (for
instance velocity vectors) to filament nodes.

See also the [VTK documentation](https://examples.vtk.org/site/VTKFileFormats/#hdf-file-formats)
for details on the VTK HDF format.

## Typical usage

```julia
using HDF5
h5open("filaments.hdf", "w") do io
    init_vtkhdf(io, fs)                   # here `fs` is a list of filaments
    write_point_data(io, "Velocity", vs)  # here `vs` contains one velocity vector per filament node
    write_field_data(io, "Time", 0.3)
    # one can add other fields here...
end
```
"""
function init_vtkhdf(
        io::HDF5.File,
        fs::AbstractVector{<:AbstractFilament},
    )
    gtop = HDF5.create_group(io, "VTKHDF")
    HDF5.attrs(gtop)["Version"] = [1, 0]

    # This attribute *must* be written as ASCII instead of UTF8, or ParaView will fail to
    # open the file.
    let s = "UnstructuredGrid"
        dtype = datatype_ascii(s)
        dspace = HDF5.dataspace(s)
        attr = HDF5.create_attribute(gtop, "Type", dtype, dspace)
        HDF5.write_attribute(attr, dtype, s)
    end

    V = eltype(eltype(fs))  # usually V == SVector{3, T}
    T = eltype(V)
    N = length(V)  # usually == 3

    num_points = sum(f -> length(nodes(f)) + 1, fs)  # the +1 is to include the endpoint
    num_cells = length(fs)

    gtop["NumberOfCells"] = [num_cells]
    gtop["NumberOfPoints"] = [num_points]
    gtop["NumberOfConnectivityIds"] = [num_points]

    dset_points = HDF5.create_dataset(gtop, "Points", T, (N, num_points))
    dset_connec = HDF5.create_dataset(gtop, "Connectivity", Int, (num_points,))
    dset_offset = HDF5.create_dataset(gtop, "Offsets", Int, (num_cells + 1,))
    dset_ctypes = HDF5.create_dataset(gtop, "Types", UInt8, (num_cells,))

    let n = 0
        dset_offset[1] = 0
        for (i, f) ∈ enumerate(fs)
            Xs = @views nodes(f)[begin:(end + 1)]
            Np = length(Xs)
            dset_points[:, (n + 1):(n + Np)] = reinterpret(reshape, T, Xs)
            dset_connec[(n + 1):(n + Np)] = n:(n + Np - 1)  # switch to zero-based indexing
            dset_ctypes[i] = UInt8(4)  # cell type, `4` corresponds to VTK_POLY_LINE
            n += Np
            dset_offset[i + 1] = n
        end
        @assert n == num_points
    end

    close(gtop)
    nothing
end

"""
    write_point_data(
        io::HDF5.File,
        name::AbstractString,
        vs::AbstractVector{<:AbstractVector},
    )

Attach data to filament nodes.

This can be used to write fields defined at filament nodes (for instance, the velocity of
each node).

The data is written to the dataset `/VTKHDF/PointData/\$name`.

For vector fields (such as velocity), the written dataset has dimensions `(3, Np)` where
`Np` is the total number of filament nodes (including endpoints).
The format is exactly the same as for the `Points` dataset as detailed in
[`init_vtkhdf`](@ref). As also explained there, the `Offsets` dataset can be used to recover
the values associated to each filament.

Note that [`init_vtkhdf`](@ref) must be called *before* using this function.
"""
function write_point_data(
        io::HDF5.File,
        name::AbstractString,
        vs::AbstractVector{<:AbstractVector},
    )
    gtop = HDF5.open_group(io, "VTKHDF")
    gdata = open_or_create_group(gtop, "PointData")

    # For now, assume `vs` contains vector data (such as velocities).
    # This can be easily generalised later.
    attrname = "Vectors"
    names = if haskey(HDF5.attrs(gdata), attrname)
        names_prev = HDF5.attrs(gdata)[attrname] :: Vector{String}
        vcat(names_prev, name) :: Vector{String}
    else
        String[name]
    end
    HDF5.attrs(gdata)[attrname] = names  # overwrite attribute if it already existed

    # Write actual data
    V = eltype(eltype(vs))  # usually V == SVector{3, T}
    T = eltype(V)
    N = length(V)  # usually == 3
    num_points = sum(v -> length(v) + 1, vs)  # the +1 is to include the endpoint
    let dset = HDF5.create_dataset(gdata, name, T, (N, num_points))
        n = 0
        for v ∈ vs
            # data = @view v[begin:(end + 1)]
            # Np = length(data)
            Np = length(v)
            dset[:, (n + 1):(n + Np)] = reinterpret(reshape, T, v)
            dset[:, n + Np + 1] = v[begin]  # assume end point == start point
            n += Np + 1
        end
        @assert n == num_points
        close(dset)
    end

    close(gdata)
    close(gtop)

    nothing
end

"""
    write_field_data(io::HDF5.File, name::AbstractString, data)

Write data as VTK field data to VTK HDF file.

In VTK, *field data* refers to data which is not directly attached to the geometry.
This is typically small datasets or simple values, such as the current time or simulation
parameters.
"""
function write_field_data(io::HDF5.File, name::AbstractString, data)
    gtop = HDF5.open_group(io, "VTKHDF")
    gdata = open_or_create_group(gtop, "FieldData")
    _write_field_data(gdata, name, data)
    close(gdata)
    close(gtop)
    nothing
end

# General arrays of data (this actually will only work with regular `Array`s I think...)
function _write_field_data(g, name, data::AbstractArray)
    g[name] = data
end

# Strings as VTK field data require special handling:
# - they must be encoded in ASCII
# - they must be variable-length strings (default in HDF5.jl for arrays of strings, but we
#   enforce it anyways just in case)
function _write_field_data(g, name, data::AbstractArray{<:AbstractString})
    dtype = HDF5.datatype(data)
    HDF5.API.h5t_set_cset(dtype.id, HDF5.API.H5T_CSET_ASCII)
    HDF5.API.h5t_set_size(dtype.id, HDF5.API.H5T_VARIABLE)
    dspace = HDF5.dataspace(data)
    dset = HDF5.create_dataset(g, name, dtype, dspace)
    HDF5.write_dataset(dset, dtype, data)
end

# Otherwise, assume scalar data (numbers, strings)
function _write_field_data(g, name, data)
    _write_field_data(g, name, [data])
end

"""
    read_filaments(io::HDF5.File, ::Type{T}, method::DiscretisationMethod) -> Vector{<:AbstractFilament}

Read filament locations from VTK HDF file.

This function loads filaments based on the datasets `/VTKHDF/Points` and `/VTKHDF/Offsets`
as written by the [`init_vtkhdf`](@ref) function.

Returns a vector of filaments. The specific type of filament (e.g.
[`ClosedSplineFilament`](@ref) or [`ClosedLocalFilament`](@ref)) depends on the chosen
`method`. See [`Filaments.init`](@ref) for possible options.

## Typical usage

```
using HDF5
h5open("filaments.hdf", "r") do io
    fs = read_filaments(io, Float64, CubicSplineMethod())  # here `fs` is a list of filaments
    vs = read_point_data(io, "Velocity")  # here `vs` contains one velocity vector per filament node
    t = read_field_data(io, "Time")
    # one can read other fields here...
end
```
"""
function read_filaments(
        io::HDF5.File,
        ::Type{T},
        method::DiscretisationMethod,
    ) where {T <: AbstractFloat}
    haskey(io, "VTKHDF") || error("expected a `VTKHDF` group at the top of the file")
    gtop = HDF5.open_group(io, "VTKHDF")

    dset_points = HDF5.open_dataset(gtop, "Points")
    dset_offset = HDF5.open_dataset(gtop, "Offsets")

    Ñ, num_points = size(dset_points) :: Dims{2}
    @assert Ñ == 3
    num_cells = length(dset_offset) - 1
    num_cells > 0 || error("expected at least one filament in the file")

    # Read points as an Array{T} even if it was written as an Array{S} with S ≠ T
    # (for instance, if data was written as S = Float64 and we want to read it as T = Float32).
    # Similar thing for offsets. This ensures type stability...
    points = Array{T}(undef, 3, num_points)
    offsets = Array{Int}(undef, num_cells + 1)

    HDF5.read_dataset(dset_points, HDF5.datatype(dset_points), points)
    HDF5.read_dataset(dset_offset, HDF5.datatype(dset_offset), offsets)

    fs = map(1:num_cells) do i
        _load_filament(T, points, offsets, method, i)
    end

    close(gtop)

    fs
end

# Load filament `i` from the list of points using the given offsets
function _load_filament(
        ::Type{T}, points::AbstractMatrix, offsets::AbstractVector, method, i::Int,
    ) where {T}
    S = eltype(points)
    Xs = reinterpret(Vec3{S}, points)
    a = offsets[i]      # points[:, a + 1] is the first point of this filament
    b = offsets[i + 1]  # points[:, b] is the filament endpoint
    num_nodes = b - a - 1         # excluding the endpoint
    Xoffset = Xs[b] - Xs[a + 1]
    f = Filaments.init(ClosedFilament{T}, num_nodes, method; offset = Xoffset)
    @inbounds for j ∈ eachindex(f)
        f[j] = Xs[a + j]
    end
    f
end
