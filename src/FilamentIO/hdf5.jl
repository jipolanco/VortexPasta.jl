export init_vtkhdf, write_point_data

using HDF5: HDF5

## ================================================================================ ##
## 1. Writing filament HDF5 files
## ================================================================================ ##

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
    init_vtkhdf(io::HDF5.File, fs::AbstractVector{<:AbstractFilament}; refinement = 1)

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

## Optional arguments

- `refinement::Int = 1`: allows to output more than 1 point for each filament segment. This
  is mostly useful for producing nice visualisations. The level of refinement is writen to
  the `/VTKHDF/RefinementLevel` dataset, which allows to read back the data skipping
  intra-segment nodes.

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
        fs::AbstractVector{<:AbstractFilament};
        refinement::Int = 1,
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

    num_points = sum(f -> refinement * length(nodes(f)) + 1, fs)  # the +1 is to include the endpoint
    num_cells = length(fs)

    gtop["NumberOfCells"] = [num_cells]
    gtop["NumberOfPoints"] = [num_points]
    gtop["NumberOfConnectivityIds"] = [num_points]
    gtop["RefinementLevel"] = refinement  # this is not a VTK attribute, it's just for our own use

    points = let
        dtype = HDF5.datatype(T) 
        dspace = HDF5.dataspace((N, num_points))
        (;
            dtype, dspace,
            dset = HDF5.create_dataset(gtop, "Points", dtype, dspace),
        )
    end

    let dset = HDF5.create_dataset(gtop, "Connectivity", Int, (num_points,))
        local data = collect(0:(num_points - 1))
        dtype = HDF5.datatype(data)
        HDF5.write_dataset(dset, dtype, data)
        close(dtype)
        close(dset)
    end

    let dset = HDF5.create_dataset(gtop, "Offsets", Int, (num_cells + 1,))
        local data = Vector{Int}(undef, num_cells + 1)
        data[1] = 0
        for (i, f) ∈ enumerate(fs)
            data[i + 1] = data[i] + refinement * length(nodes(f)) + 1
        end
        dtype = HDF5.datatype(data)
        HDF5.write_dataset(dset, dtype, data)
        close(dtype)
        close(dset)
    end

    # Write cell types, `4` corresponds to VTK_POLY_LINE
    let dset = HDF5.create_dataset(gtop, "Types", UInt8, (num_cells,))
        local data = fill(UInt8(4), num_cells)
        HDF5.write_dataset(dset, HDF5.datatype(data), data)
        close(dset)
    end

    n = 0
    for f ∈ fs
        if refinement == 1
            Xs = nodes(f)
            Np = length(Xs) + 1  # include the endpoint
        else
            # Reuse Makie recipe code.
            Xs = Filaments._refine_filament(f, refinement)
            @assert typeof(Xs) === typeof(nodes(f))
            @assert Xs[end] ≈ f[end + 1]  # already includes the endpoint
            Np = length(Xs)
            @assert Np == refinement * length(nodes(f)) + 1
        end

        let (; dset, dspace, dtype,) = points
            memspace = HDF5.dataspace((N, Np))
            memtype = dtype
            HDF5.select_hyperslab!(dspace, (1:N, (n + 1):(n + Np)))
            HDF5.API.h5d_write(dset.id, memtype.id, memspace.id, dspace.id, dset.xfer, Xs)
            close(memspace)
        end

        n += Np
    end
    @assert n == num_points

    map(close, points)  # close dataset + dataspace + datatype

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

    refinement = HDF5.read_dataset(gtop, "RefinementLevel") :: Int

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
    num_points = sum(v -> refinement * length(v) + 1, vs)  # the +1 is to include the endpoint
    dtype = HDF5.datatype(T)
    dspace = HDF5.dataspace((N, num_points))
    memtype = dtype
    memspace_vec = HDF5.dataspace((N,))  # this is for writing one vector at a time

    dset = HDF5.create_dataset(gdata, name, dtype, dspace)
    n = 0
    for v ∈ vs
        Np = refinement * length(v)
        if refinement == 1
            # Write all data associated to a single filament at once
            memspace_filament = HDF5.dataspace((N, Np))
            HDF5.select_hyperslab!(dspace, (1:N, (n + 1):(n + Np)))
            HDF5.API.h5d_write(dset.id, memtype.id, memspace_filament.id, dspace.id, dset.xfer, v)
            close(memspace_filament)
            n += Np
        else
            # For now just use linear interpolation
            for i ∈ eachindex(v), m ∈ 1:refinement
                α = (m - 1) / refinement
                u = (1 - α) * v[i] + α * v[i + 1]  # assumes v is padded!!
                n += 1
                # Write vectors one by one
                HDF5.select_hyperslab!(dspace, (1:N, n))
                HDF5.API.h5d_write(dset.id, memtype.id, memspace_vec.id, dspace.id, dset.xfer, u)
            end
        end
        n += 1
        HDF5.select_hyperslab!(dspace, (1:N, n))
        HDF5.API.h5d_write(dset.id, memtype.id, memspace_vec.id, dspace.id, dset.xfer, v[begin])
    end
    @assert n == num_points
    close(dset)

    close(dspace)
    close(dtype)
    close(memspace_vec)
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

Note that scalar data (such as time) is always written as a single-element vector, since
otherwise it cannot be parsed by VTK.
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
# - they must be variable-length strings
function _write_field_data(g, name, data::AbstractArray{<:AbstractString})
    # This is adapted from the code for HDF5.datatype(data), but making sure we use ASCII.
    S = eltype(data)
    type_id = HDF5.API.h5t_copy(HDF5.hdf5_type_id(S))
    HDF5.API.h5t_set_size(type_id, HDF5.API.H5T_VARIABLE)
    HDF5.API.h5t_set_cset(type_id, HDF5.API.H5T_CSET_ASCII)
    dtype = HDF5.Datatype(type_id)
    dspace = HDF5.dataspace(data)
    dset = HDF5.create_dataset(g, name, dtype, dspace)
    HDF5.write_dataset(dset, dtype, data)
end

# Otherwise, assume scalar data (numbers, strings)
function _write_field_data(g, name, data)
    _write_field_data(g, name, [data])
end

## ================================================================================ ##
## 2. Reading filament HDF5 files
## ================================================================================ ##

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
    vs = read_point_data(io, "Velocity", fs)        # here `vs` contains one velocity vector per filament node
    t = only(read_field_data(io, "Time", Float64))  # note: field data is always written as an array
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

    refinement = HDF5.read_dataset(gtop, "RefinementLevel") :: Int

    dset_points = HDF5.open_dataset(gtop, "Points")
    dset_offset = HDF5.open_dataset(gtop, "Offsets")

    Ñ, num_points_refined = size(dset_points) :: Dims{2}
    @assert Ñ == 3
    num_cells = length(dset_offset) - 1
    num_cells > 0 || error("expected at least one filament in the file")

    if refinement == 1
        num_points = num_points_refined
    else
        # We subtract endpoint from each cell before dividing by refinement level
        num_points, _remainder = divrem(num_points_refined - num_cells, refinement)
        @assert _remainder == 0  # the division is exact
        num_points += num_cells  # reinclude the endpoints in final count
    end

    # Read points as an Array{T} even if it was written as an Array{S} with S ≠ T
    # (for instance, if data was written as S = Float64 and we want to read it as T = Float32).
    # Similar thing for offsets. This ensures type stability...
    points = Array{T}(undef, 3, num_points)
    offsets = Array{Int}(undef, num_cells + 1)

    read_dataset!(offsets, dset_offset)

    if refinement == 1
        read_dataset!(points, dset_points)
    else
        _recompute_offsets!(offsets, refinement)
        @assert last(offsets) == num_points
        dspace = HDF5.dataspace(dset_points)
        n = 0
        for i ∈ 1:num_cells
            a = offsets[i] + 1
            b = offsets[i + 1]
            out = @view points[:, a:b]
            Np = size(out, 2)
            inds_file = range(n + 1; length = Np, step = refinement)
            HDF5.select_hyperslab!(dspace, (1:Ñ, inds_file))
            read_dataset!(out, dset_points, dspace)
            n = last(inds_file)
        end
    end

    fs = map(1:num_cells) do i
        _load_filament(T, points, offsets, method, i)
    end

    close(gtop)

    fs
end

function _load_with_refinement!()
end

function _recompute_offsets!(offsets, refinement)
    Base.require_one_based_indexing(offsets)
    num_cells = length(offsets) - 1
    off_i = offsets[begin]  # offset of previous cell before recomputation
    for i ∈ 1:num_cells
        Np_refined = offsets[i + 1] - off_i
        Np, _remainder = divrem(Np_refined - 1, refinement)
        @assert _remainder == 0
        Np += 1  # reinclude the endpoint in final count
        off_i = offsets[i + 1]
        offsets[i + 1] = offsets[i] + Np
    end
    offsets
end

# Load filament `i` from the list of points using the given offsets
function _load_filament(
        ::Type{T}, points::AbstractMatrix, offsets::AbstractVector, method, i::Int,
    ) where {T}
    S = eltype(points)
    Xs = reinterpret(Vec3{S}, points)
    a = offsets[i]      # points[:, a + 1] is the first point of this filament
    b = offsets[i + 1]  # points[:, b] is the filament endpoint
    num_nodes = b - a - 1  # not counting the endpoint
    Xoffset = Xs[b] - Xs[a + 1]
    f = Filaments.init(ClosedFilament{T}, num_nodes, method; offset = Xoffset)
    @inbounds for j ∈ eachindex(f)
        f[j] = Xs[a + j]
    end
    f
end

"""
    read_point_data(
        io::HDF5.File, name::AbstractString,
        fs::AbstractVector{<:AbstractFilament},
        [V = eltype(eltype(fs))],  # default is Vec3{T}
    ) -> AbstractVector{<:AbstractVector}

Load point data from VTK HDF file.

Reads data written by [`write_point_data`](@ref).
The output is a vector of vectors, where each subvector holds the data associated to a
single filament.

See [`read_filaments`](@ref) for a typical usage example.
"""
function read_point_data(
        io::HDF5.File, name::AbstractString,
        fs::AbstractVector{<:AbstractFilament},
        ::Type{V} = eltype(eltype(fs)),  # Vec3{T} by default
    ) where {V}
    vs = map(f -> similar(nodes(f), V), fs)  # one PaddedVector for each filament
    read_point_data!(io, vs, name)
end

"""
    read_point_data!(
        io::HDF5.File,
        vs::AbstractVector{<:AbstractVector},
        name::AbstractString,
    )

Load point data from VTK HDF file.

Data is read onto the vector of vectors `vs`.

See also [`read_point_data`](@ref) and [`write_point_data`](@ref).
"""
function read_point_data!(
        io::HDF5.File,
        vs::AbstractVector{<:AbstractVector},
        name::AbstractString,
    )
    # TODO generalise for scalar data
    V = eltype(eltype(vs))  # usually V == SVector{3, T}
    N = length(V)  # usually == 3
    T = eltype(V)

    gtop = HDF5.open_group(io, "VTKHDF")
    refinement = HDF5.read_dataset(gtop, "RefinementLevel") :: Int

    gdata = HDF5.open_group(gtop, "PointData")
    dset = HDF5.open_dataset(gdata, name)
    Ñ, num_points = size(dset) :: Dims{2}
    @assert N == Ñ

    n = 0
    dspace = HDF5.dataspace(dset)
    for v ∈ vs
        Np = length(v)  # number of filament nodes (*excluding* the endpoint)
        vdata = reinterpret(reshape, T, v)
        if refinement == 1
            HDF5.select_hyperslab!(dspace, (1:N, (n + 1):(n + Np)))  # read part of the dataset
            read_dataset!(vdata, dset, dspace)
            n += Np + 1  # the + 1 is because the endpoint is included in the file (but ignored here)
        else
            inds_file = range(n + 1; length = Np, step = refinement)
            HDF5.select_hyperslab!(dspace, (1:N, inds_file))
            read_dataset!(vdata, dset, dspace)
            n = last(inds_file) + refinement
        end
    end
    @assert n == num_points

    close(dset)
    close(gdata)
    close(gtop)

    vs
end

function read_dataset!(
        out::AbstractArray,
        dset::HDF5.Dataset,
        dspace = HDF5.dataspace(dset);  # this can also be a hyperslab
        memspace = HDF5.dataspace(out),
    )
    memtype = HDF5.datatype(out)
    HDF5.API.h5d_read(dset, memtype, memspace, dspace, dset.xfer, out)
    out
end

"""
    read_field_data(
        io::HDF5.File, name::AbstractString, ::Type{T},
    ) -> Vector{T}

Read field data from VTK HDF file.

Returns a *vector* of elements of type `T`, because, in VTK, field data is always
represented as an array (even for single values such as time).

See [`write_field_data`](@ref) for more details.
"""
function read_field_data(
        io::HDF5.File, name::AbstractString, ::Type{T},
    ) where {T}
    dset = HDF5.open_dataset(io, "/VTKHDF/FieldData/$name")
    N = length(dset)
    out = Vector{T}(undef, N)
    memtype = HDF5.datatype(out)
    memspace = HDF5.dataspace(out)
    dspace = HDF5.dataspace(dset)
    HDF5.API.h5d_read(dset, memtype, memspace, dspace, dset.xfer, out)
    close(dset)
    out
end

# Specific case of strings
function read_field_data(
        io::HDF5.File, name::AbstractString, ::Type{<:AbstractString},
    )
    dset = HDF5.open_dataset(io, "/VTKHDF/FieldData/$name")
    read(dset) :: Vector{String}
end
