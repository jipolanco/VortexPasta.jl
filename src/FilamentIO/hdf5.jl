export init_vtkhdf, write_point_data

using HDF5: HDF5

# This is to make sure we write the dataset type ("UnstructuredGrid") as ASCII.
# This is required by VTK/ParaView.
function datatype_ascii(s::AbstractString)
    dtype = HDF5.datatype(s)
    HDF5.API.h5t_set_cset(dtype.id, HDF5.API.H5T_CSET_ASCII)
    dtype
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
  Once again, note that number of nodes of filament `i` *excluding the endpoint* is
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

    gname = "PointData"
    gdata = if haskey(gtop, gname)
        HDF5.open_group(gtop, gname)
    else
        HDF5.create_group(gtop, gname)
    end

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
