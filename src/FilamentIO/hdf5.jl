export write_vtkhdf, read_vtkhdf, PointData, FieldData

using HDF5: HDF5

## ================================================================================ ##
## 1. Writing filament HDF5 files
## ================================================================================ ##

struct VTKHDFFile{
        Group     <: HDF5.Group,
        Filaments <: AbstractVector{<:AbstractFilament},
    }
    gtop       :: Group  # "/VTKHDF" group
    fs         :: Filaments
    refinement :: Int
end

struct PointData end
struct FieldData end

"""
    write_vtkhdf(
        [func::Function],
        filename::AbstractString,
        fs::AbstractVector{<:AbstractFilament};
        refinement = 1,
    )

Write new VTK HDF file containing a list of filaments.

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

See also the [VTK documentation](https://examples.vtk.org/site/VTKFileFormats/#hdf-file-formats)
for details on the VTK HDF format.

## Attaching extra data

The optional `func` argument can be used to attach other data (such as velocity vectors or
the current time) to the generated file. This is most conveniently done using the `do` block
syntax. See further below for some examples.

## Optional arguments

- `refinement::Int = 1`: allows to output more than 1 point for each filament segment. This
  is mostly useful for producing nice visualisations. The level of refinement is writen to
  the `/VTKHDF/RefinementLevel` dataset, which allows to read back the data skipping
  intra-segment nodes.

## Typical usage

```julia
# Note: the file extension is arbitrary, but ParaView prefers ".hdf" if one wants to use the
# files for visualisation.
write_vtkhdf("filaments.hdf", fs; refinement = 2) do io
    io["Velocity"] = vs  # adds velocity as VTK point data, assuming vs is a VectorOfVectors
    io["Time"] = 0.3     # adds time as VTK field data, since it's a scalar
    # one can add other fields here...
end
```
"""
function write_vtkhdf(
        func::Func,
        filename::AbstractString,
        fs::AbstractVector{<:AbstractFilament};
        refinement = 1,
    ) where {Func <: Function}
    HDF5.h5open(filename, "w") do io
        writer = init_vtkhdf(io, fs; refinement)
        func(writer)
    end
    nothing
end

write_vtkhdf(filename::AbstractString, args...; kws...) =
    write_vtkhdf(Returns(nothing), filename, args...; kws...)

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

    points = let dtype = HDF5.datatype(T)
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

    VTKHDFFile(gtop, fs, refinement)
end

vtk_attribute_type(::Type{<:Number}) = "Scalars"          # usually <:AbstractFloat
vtk_attribute_type(::Type{<:AbstractVector}) = "Vectors"  # usually Vec3

"""
    Base.setindex!(io::VTKHDFFile, vs::AbstractVector{<:AbstractVector}, name::AbstractString)

Attach data to filament nodes.

One generally wants to use the syntax `io[name] = vs` which calls this function.

This can be used to write fields defined at filament nodes (for instance, the velocity of
each node).

The data is written to the dataset `/VTKHDF/PointData/\$name`.

For vector fields (such as velocity), the written dataset has dimensions `(3, Np)` where
`Np` is the total number of filament nodes (including endpoints).
The format is exactly the same as for the `Points` dataset as detailed in
[`write_vtkhdf`](@ref). As also explained there, the `Offsets` dataset can be used to recover
the values associated to each filament.
"""
function Base.setindex!(
        writer::VTKHDFFile,
        vs::AbstractVector{<:AbstractVector},
        name::AbstractString,
    )
    (; gtop, refinement,) = writer
    gdata = open_or_create_group(gtop, "PointData")
    T = eltype(eltype(vs))  # usually <:Real or Vec3{<:Real}
    attrname = vtk_attribute_type(T)  # "Scalars", "Vectors", ...
    names = if haskey(HDF5.attrs(gdata), attrname)
        names_prev = HDF5.attrs(gdata)[attrname] :: Vector{String}
        vcat(names_prev, name) :: Vector{String}
    else
        String[name]
    end
    HDF5.attrs(gdata)[attrname] = names  # overwrite attribute if it already existed
    _write_data_on_filaments(T, gdata, vs, name, refinement)
    close(gdata)
    nothing
end

# Write data on filament nodes -- case of scalar data (e.g. local curvature magnitude)
function _write_data_on_filaments(
        ::Type{T}, gdata, vs::AbstractVector{<:AbstractVector}, name, refinement,
    ) where {T <: Number}
    @assert T === eltype(eltype(vs))
    num_points = sum(v -> refinement * length(v) + 1, vs)  # the +1 is to include the endpoint
    dtype = HDF5.datatype(T)
    dspace = HDF5.dataspace((num_points,))
    memtype = dtype
    memspace_node = HDF5.dataspace(T)  # this is for writing a scalar on a single point

    dset = HDF5.create_dataset(gdata, name, dtype, dspace)
    n = 0
    for v ∈ vs
        Np = refinement * length(v)
        if refinement == 1
            # Write all data associated to a single filament at once
            memspace_filament = HDF5.dataspace((Np,))
            HDF5.select_hyperslab!(dspace, ((n + 1):(n + Np),))
            HDF5.API.h5d_write(dset.id, memtype.id, memspace_filament.id, dspace.id, dset.xfer, v)
            close(memspace_filament)
            n += Np
        else
            # For now just use linear interpolation
            for i ∈ eachindex(v), m ∈ 1:refinement
                u = _interpolate_on_segment(v, i, m, refinement)
                n += 1
                # Write values one by one
                HDF5.select_hyperslab!(dspace, (n,))
                HDF5.API.h5d_write(dset.id, memtype.id, memspace_node.id, dspace.id, dset.xfer, Ref(u))
            end
        end
        n += 1
        HDF5.select_hyperslab!(dspace, (n,))
        HDF5.API.h5d_write(dset.id, memtype.id, memspace_node.id, dspace.id, dset.xfer, v)
    end
    @assert n == num_points
    close(dset)

    close(dspace)
    close(dtype)
    close(memspace_node)

    nothing
end

# Write data on filament nodes -- case of vector data (e.g. velocity)
function _write_data_on_filaments(
        ::Type{V}, gdata, vs::AbstractVector{<:AbstractVector}, name, refinement,
    ) where {V <: AbstractVector}
    @assert V === eltype(eltype(vs))
    T = eltype(V)
    @assert T <: Number
    N = length(V)  # usually == 3
    num_points = sum(v -> refinement * length(v) + 1, vs)  # the +1 is to include the endpoint
    dtype = HDF5.datatype(T)
    dspace = HDF5.dataspace((N, num_points))
    memtype = dtype
    memspace_node = HDF5.dataspace((N,))  # this is for writing a vector on a single point

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
                u = _interpolate_on_segment(v, i, m, refinement)
                n += 1
                # Write vectors one by one
                HDF5.select_hyperslab!(dspace, (1:N, n))
                HDF5.API.h5d_write(dset.id, memtype.id, memspace_node.id, dspace.id, dset.xfer, u)
            end
        end
        n += 1
        HDF5.select_hyperslab!(dspace, (1:N, n))
        HDF5.API.h5d_write(dset.id, memtype.id, memspace_node.id, dspace.id, dset.xfer, v[begin])
    end
    @assert n == num_points
    close(dset)

    close(dspace)
    close(dtype)
    close(memspace_node)

    nothing
end

# Interpolate value on segment. For now we simply use linear interpolation (it's just for
# visualisation anyways).
function _interpolate_on_segment(v, i, m, M)
    α = (m - 1) / M
    (1 - α) * v[i] + α * v[i + 1]  # assumes v is padded!!
end

"""
    Base.setindex!(io::VTKHDFFile, data, name::AbstractString)

Write data as VTK field data to VTK HDF file.

In VTK, *field data* refers to data which is not directly attached to the geometry.
This is typically small datasets or simple values, such as the current time or simulation
parameters.

Note that scalar data (such as time) is always written as a single-element vector, since
otherwise it cannot be parsed by VTK.

This function interprets everything that is not a vector of vectors as field data.
"""
function Base.setindex!(writer::VTKHDFFile, data, name::AbstractString)
    (; gtop,) = writer
    gdata = open_or_create_group(gtop, "FieldData")
    _write_field_data(gdata, name, data)
    close(gdata)
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
    read_vtkhdf(
        [func::Function],
        filename::AbstractString,
        ::Type{T},
        method::DiscretisationMethod,
    ) where {T}

Read filament locations from VTK HDF file.

This function loads filaments based on the datasets `/VTKHDF/Points` and `/VTKHDF/Offsets`
as written by the [`write_vtkhdf`](@ref) function.

Returns a vector of filaments. The specific type of filament (e.g.
[`ClosedSplineFilament`](@ref) or [`ClosedLocalFilament`](@ref)) depends on the chosen
`method`. See [`Filaments.init`](@ref) for possible options.

One can also read other datasets using `read` and `read!`, as shown and explained below.

## Typical usage

```julia
local vs, t  # make sure these variables still exist after the `do` block

# The returned `fs` is a list of filaments.
fs = read_vtkhdf("filaments.hdf", Float64, CubicSplineMethod()) do io
    vs = read(io, "Velocity", PointData())            # here `vs` contains one velocity vector per filament node
    t = only(read(io, "Time", FieldData(), Float64))  # note: field data is always written as an array
    # one can read other fields here...
end
```

The available reading functions are:

- `read(io, name::AbstractString, ::PointData)` for reading point data (e.g. a velocity field);

- `read(io, name::AbstractString, ::FieldData, ::Type{T})` for reading field data (i.e. data
  not attached to filament nodes, such as the current time);

- `read!(io, vs::AbstractVector{<:AbstractVector}, name::AbstractString)` for reading point
  data onto a preallocated vector of vectors.
"""
function read_vtkhdf(
        func::Func, filename::AbstractString, ::Type{T},
        method::DiscretisationMethod,
    ) where {Func <: Function, T}
    HDF5.h5open(filename, "r") do io
        reader = read_filaments(io, T, method)
        func(reader)
        reader.fs
    end
end

read_vtkhdf(filename::AbstractString, args...; kws...) =
    read_vtkhdf(Returns(nothing), filename, args...; kws...)

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

    VTKHDFFile(gtop, fs, refinement)
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
    update_coefficients!(f)  # compute interpolations and derivative coefficients
    f
end

function Base.read(
        reader::VTKHDFFile, name::AbstractString, ::PointData,
        ::Type{V} = eltype(eltype(reader.fs)),  # Vec3{T} by default
    ) where {V}
    (; fs,) = reader
    vs = map(f -> similar(nodes(f), V), fs)  # one PaddedVector for each filament
    read!(reader, vs, name)
end

function Base.read!(
        reader::VTKHDFFile,
        vs::AbstractVector{<:AbstractVector},
        name::AbstractString,
    )
    # TODO generalise for scalar data
    V = eltype(eltype(vs))  # usually V == SVector{3, T}
    N = length(V)  # usually == 3
    T = eltype(V)

    (; gtop, refinement,) = reader

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
        if v isa PaddedVector
            pad_periodic!(v)
        end
    end
    @assert n == num_points

    close(dset)
    close(gdata)

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

function Base.read(
        reader::VTKHDFFile, name::AbstractString, ::FieldData, ::Type{T},
    ) where {T}
    dset = HDF5.open_dataset(reader.gtop, "FieldData/$name")
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
function Base.read(
        reader::VTKHDFFile, name::AbstractString, ::FieldData, ::Type{<:AbstractString},
    )
    dset = HDF5.open_dataset(reader.gtop, "FieldData/$name")
    read(dset) :: Vector{String}
end
