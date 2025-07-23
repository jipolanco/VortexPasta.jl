export write_vtkhdf, read_vtkhdf, PointData, FieldData

using HDF5: HDF5
using Bumper: Bumper, @no_escape, @alloc

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

root(io::VTKHDFFile) = HDF5.root(io.gtop)::HDF5.Group  # returns the "/" group

const FilamentId = Int32
const FilamentPeriodicOffset = Int32

struct PointData end
struct FieldData end

include("json_vtk_series.jl")
include("write_points.jl")
include("write_parametrisation.jl")
include("write_data.jl")

# Here "padded" means that Xs[end] corresponds the endpoint (which is equal to the start
# point for closed curves).
function write_padded_data_to_hdf5(points, Xs, n, ldims, Np)
    (; dset, dspace, dtype,) = points
    memspace = HDF5.dataspace((ldims..., Np))
    memtype = dtype
    lranges = map(N -> 1:N, ldims)
    HDF5.select_hyperslab!(dspace, (lranges..., (n + 1):(n + Np)))
    HDF5.API.h5d_write(dset.id, memtype.id, memspace.id, dspace.id, dset.xfer, Xs)
    close(memspace)
    nothing
end

"""
    write_vtkhdf(
        [f::Function],
        filename::AbstractString,
        fs::AbstractVector{<:AbstractFilament};
        refinement = 1,
        dataset_type = :PolyData,
        parametrisation = true,
        periods = nothing,
    )

Write new VTK HDF file containing a list of filaments.

A VTK HDF file is an HDF5 file organised in such a way that it can be readily visualised in
tools such as ParaView.
Data in the file can be readily accessed using HDF5 tools and libraries.

Following the VTK HDF specification, this function creates a "VTKHDF" group on top of the
HDF5 file. Then, it creates the datasets allowing to describe the filaments as an
unstructured grid.

See also the [VTK documentation](https://docs.vtk.org/en/latest/design_documents/VTKFileFormats.html#vtkhdf-file-format)
for details on the VTK HDF format.

# Extended help

## Relevant datasets

Some relevant datasets which are written are:

- `/VTKHDF/Points`: contains the coordinates of all filament nodes. Points are represented as
  an array of dimensions `(3, Np)` where `Np` is the total number of nodes. Note that the
  points *include the endpoint*, meaning that for a closed filament the initial coordinate
  appears twice in the file. This is done to disambiguate between *closed* and *infinite but
  unclosed* filaments (see [`end_to_end_offset`](@ref)).

- `/VTKHDF/CellData/FilamentIds`: contains the filament id associated to each VTK cell. It is
  a dataset of length `Nc`, where `Nc` is the number of cells.
  In our case a cell is a spatial curve made of discrete points.
  The ids are `$FilamentId` values from `1` to `Nf`, where `Nf` is the number of filaments.
  The values are always sorted increasingly (`FilamentIds[i + 1] ≥ FilamentIds[i]` for all `i`).
  Note that, in general, one cell corresponds to one filament (so that `Nc == Nf`), but this
  may not be the case when the `periods` argument is used and a filament is broken onto
  multiple curves (see **Periodic wrapping of filaments** below for details).

## Attaching extra data

The optional `f` argument can be used to attach other data (such as velocity vectors or
the current time) to the generated file. This is most conveniently done using the `do` block
syntax. See further below for some examples.

## Optional keyword arguments

- `refinement::Int = 1`: allows to output more than 1 point for each filament segment. This
  is mostly useful for producing nice visualisations. The level of refinement is written to
  the `/VTKHDF/RefinementLevel` dataset, which allows to read back the data skipping
  intra-segment nodes.

- `parametrisation::Bool = true`: if `true` (default), write curve parametrisation values.
  This allows `read_vtkhdf` to reconstruct the exact same curve that was written (assuming
  the same discretisation method is used, e.g. cubic splines), even if the curve used some
  non-default curve parametrisation.
  Values are written to the `/VTKHDF/PointData/Parametrisation` dataset.

- `dataset_type::Symbol`: can be either `:PolyData` (default) or `:UnstructuredGrid`.
  There's usually no reason to change this.

## Periodic wrapping of filaments

When using periodic boundary conditions, one can use the optional `periods` argument to
periodically wrap the filaments. In this case, one should pass a tuple `periods = (Lx, Ly,
Lz)` with the period in each direction (one can pass `nothing` if there are non-periodic
directions).

In this case, filaments outside the main periodic box will be translated to fit in the
periodic box. Moreover, if a filament locally goes out of the periodic box, it will be
broken onto multiple curves so that they all fit within the domain. One can then look at the
`/VTKHDF/CellData/FilamentIds` dataset to determine which curves belong to the same
filament.

## Typical usage

```julia
# Note: the file extension is arbitrary, but ParaView prefers ".vtkhdf" (or ".hdf") if one
# wants to use the files for visualisation.
write_vtkhdf("filaments.vtkhdf", fs; refinement = 2, periods = (2π, 2π, 2π)) do io
    io["Velocity"] = vs  # adds velocity as VTK point data, assuming vs is a VectorOfVectors
    io["Curvatures"] = CurvatureVector()  # convenient syntax for writing geometric quantities along filaments
    io["Time"] = 0.3     # adds time as VTK field data, since it's a scalar
    # one can add other fields here...
end
```
"""
function write_vtkhdf(
        func::Func,
        filename::AbstractString,
        fs::AbstractVector{<:AbstractFilament};
        periods = nothing,
        kwargs...
    ) where {Func <: Function}
    if periods isa Tuple
        periods_ = periods
    else
        @assert periods isa Union{Nothing, Real}
        periods_ = (periods, periods, periods)  # same period in all directions (works with `nothing` too)
    end
    HDF5.h5open(filename, "w") do io
        writer = init_vtkhdf(io, fs; periods = periods_, kwargs...)
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

# Some internal (undocumented) arguments:
#
# - _allow_breaking_filaments: can be set to false to generate files as in old VortexPasta
#   versions (< 0.18). Used in tests only, to check that we're still able to read old VTKHDF
#   files.
#
function init_vtkhdf(
        io::HDF5.File,
        fs::AbstractVector{<:AbstractFilament};
        refinement::Int = 1,
        parametrisation = true,
        dataset_type::Symbol = :PolyData,  # either :UnstructuredGrid (old default) or :PolyData
        periods::Tuple = (nothing, nothing, nothing),
        _allow_breaking_filaments = true,
    )
    gtop = HDF5.create_group(io, "VTKHDF")
    HDF5.attrs(gtop)["Version"] = [2, 1]

    with_periods = any(!isnothing, periods)
    if with_periods && !_allow_breaking_filaments
        throw(ArgumentError("the `periods` argument requires `_allow_breaking_filaments = true`"))
    end

    # This is the first VortexPasta version needed to read back the files.
    HDF5.attrs(gtop)["VortexPasta minimum version"] = if with_periods
        [0, 18, 0]
    else
        [0, 14, 0]
    end

    # This attribute *must* be written as ASCII instead of UTF8, or ParaView will fail to
    # open the file.
    let s = string(dataset_type)
        dtype = datatype_ascii(s)
        dspace = HDF5.dataspace(s)
        attr = HDF5.create_attribute(gtop, "Type", dtype, dspace)
        HDF5.write_attribute(attr, dtype, s)
    end

    V = eltype(eltype(fs))  # usually V == SVector{3, T}
    T = eltype(V)
    N = length(V)  # usually == 3

    # Add empty categories
    if dataset_type == :PolyData
        for type ∈ ("Vertices", "Polygons", "Strips")
            gempty = HDF5.create_group(gtop, type)
            gempty["NumberOfCells"] = [0]
            gempty["NumberOfConnectivityIds"] = [0]
            gempty["Connectivity"] = Int[]
            gempty["Offsets"] = [0]
            close(gempty)
        end
    end

    num_points = sum(f -> refinement * length(nodes(f)) + 1, fs)  # the +1 is to include the endpoint

    g_pointdata = HDF5.create_group(gtop, "PointData")
    g_celldata = HDF5.create_group(gtop, "CellData")

    if dataset_type == :PolyData
        gcells = HDF5.create_group(gtop, "Lines")  # write cell info to /VTKHDF/Lines
    elseif dataset_type == :UnstructuredGrid
        gcells = gtop  # write cell info directly to /VTKHDF
    else
        throw(ArgumentError("dataset_type must be either :PolyData (default) or :UnstructuredGrid"))
    end

    gtop["NumberOfPoints"] = [num_points]

    buf = Bumper.default_buffer()
    @no_escape buf let
        # These are not VTK attributes, it's just for our own use
        gtop["RefinementLevel"] = refinement

        if _allow_breaking_filaments
            gtop["Periods"] = let
                local Ls = @alloc(T, length(periods))
                for (i, L) ∈ pairs(periods)
                    # We write a zero in "non-periodic" dimensions.
                    Ls[i] = L === nothing ? zero(T) : T(L)
                end
                Ls
            end
        end

        num_filaments = length(fs)
        gtop["NumberOfFilaments"] = [num_filaments]

        # Integer periodic offset applied to the first point of the filament.
        # This is just to be able to reconstruct the filaments with their original periodic
        # offset.
        offsets = @alloc(FilamentPeriodicOffset, N, num_filaments)

        # Points to be written.
        Xs_all = @alloc(V, num_points)

        # This determines whether the endpoint f[end + 1] of a filament `f` should be
        # included in the cell connectivity arrays. This is only useful for visualisations,
        # to avoid ugly "jumps" between points f[end] and f[end + 1] of a filament if they
        # are "far" after periodic wrapping.
        include_endpoint_in_connectivity = @alloc(Bool, num_filaments)

        (; filament_ids, cell_offsets,) = set_points!(
            Xs_all, offsets, include_endpoint_in_connectivity, fs, refinement, periods,
        )

        num_cells = length(filament_ids)
        @assert length(cell_offsets) == num_cells + 1

        # Write full Points dataset.
        points_to_vtkhdf_dataset(gtop, Xs_all)

        gcells["Offsets"] = cell_offsets

        if _allow_breaking_filaments
            g_celldata["FilamentIds"] = filament_ids
            gtop["FilamentPeriodicOffsets"] = offsets
        end

        gcells["NumberOfCells"] = [num_cells]

        # Number of points included in the connectivity vector.
        num_points_connectivity =
            num_points - (num_filaments - sum(include_endpoint_in_connectivity))
        gcells["NumberOfConnectivityIds"] = [num_points_connectivity]

        gcells["Connectivity"] = let
            connectivity = @alloc(Int, num_points_connectivity)
            local n = 0
            local m = -1  # VTK uses zero-based indexing (we start at m = 0)
            for (i, f) ∈ enumerate(fs)
                Np = length(nodes(f)) * refinement
                Np_cells = Np + include_endpoint_in_connectivity[i]
                inds = (n + 1):(n + Np_cells)
                connect = (m + 1):(m + Np_cells)
                @views copyto!(connectivity[inds], connect)
                n += Np_cells
                m += Np + 1  # always includes endpoint
            end
            @assert n == num_points_connectivity
            @assert m == num_points - 1
            connectivity
        end

        if dataset_type == :UnstructuredGrid
            # Write cell types, `4` corresponds to VTK_POLY_LINE
            let dset = HDF5.create_dataset(gtop, "Types", UInt8, (num_cells,))
                local data = @alloc(UInt8, num_cells)
                fill!(data, UInt8(4))
                HDF5.write_dataset(dset, HDF5.datatype(data), data)
                close(dset)
            end
        end

        if parametrisation
            info_param = let dtype = HDF5.datatype(T)
                dspace = HDF5.dataspace((num_points,))
                (;
                    dtype, dspace,
                    dset = HDF5.create_dataset(g_pointdata, "Parametrisation", dtype, dspace),
                )
            end
            let n = write_parametrisation(info_param, fs, refinement)
                @assert n == num_points
            end
            foreach(close, info_param)  # close dataset + dataspace + datatype
        end
    end  # @no_escape

    close(g_celldata)
    close(g_pointdata)

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

# Writing geometric quantities

It is also possible to write geometric quantities such as unit tangents or curvature vectors
along filaments (see [`GeometricQuantity`](@ref) for a list).
For this, one can pass the wanted quantity as the `vs` argument.

For example, one can simply do

    io["Curvatures"] = CurvatureVector()

to write curvature vectors along each filament.
"""
function Base.setindex!(
        writer::VTKHDFFile,
        vs::AbstractVector{<:AbstractVector},
        name::AbstractString,
    )
    write_point_data(writer, vs, name)
end

function Base.setindex!(writer::VTKHDFFile, q::GeometricQuantity, name::AbstractString)
    (; fs,) = writer
    buf = Bumper.default_buffer()
    @no_escape buf begin
        vs = map(fs) do f
            v1 = f[begin, q]
            T = typeof(v1)
            v = @alloc(T, length(f))
            v[begin] = v1
            for i ∈ eachindex(v)[2:end]
                v[i] = f[i, q]
            end
            v
        end
        write_point_data(writer, vs, name)
    end
    nothing
end

function write_point_data(writer, vs, name)
    (; gtop,) = writer
    gdata = HDF5.open_group(gtop, "PointData")
    T = eltype(eltype(vs))  # usually <:Real or Vec3{<:Real}
    attrname = vtk_attribute_type(T)  # "Scalars", "Vectors", ...
    names = if haskey(HDF5.attrs(gdata), attrname)
        names_prev = HDF5.attrs(gdata)[attrname] :: Vector{String}
        vcat(names_prev, name) :: Vector{String}
    else
        String[name]
    end
    HDF5.attrs(gdata)[attrname] = names  # overwrite attribute if it already existed
    write_data_on_filaments(T, writer, gdata, vs, name)
    close(gdata)
    nothing
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
    read_vtkhdf([f::Function], filename, T, method::DiscretisationMethod)

Read filament locations from VTK HDF file.

This function loads filaments based on the datasets `/VTKHDF/Points` and `/VTKHDF/Offsets`
as written by the [`write_vtkhdf`](@ref) function.

Returns a vector of filaments with precision `T` (`Float32` or `Float64`).
Each filament is discretised according to the chosen `method`.
See [`Filaments.init`](@ref) for possible options.

One can also read other datasets using `read` and `read!`, as shown and explained below.

# Extended help

## Typical usage

```julia
local vs, t  # make sure these variables still exist after the `do` block

# The returned `fs` is a list of filaments.
fs = read_vtkhdf("filaments.vtkhdf", Float64, CubicSplineMethod()) do io
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

## Accessing filament data in `do` block

When using the `do` block syntax as in the above example, one may want to have access to the
filament locations `fs` from within the `do` block. For instance, this can be useful if one
has a preallocated vector of velocities `vs` which needs to be resized to match the number
of filaments and filament nodes, before reading values using `read!`.

In fact, `fs` can be easily obtained from the `io` object:

```julia
read_vtkhdf("filaments.vtkhdf", Float64, CubicSplineMethod()) do io
    fs = io.fs  # this is the vector of filaments
    # Say we want to resize an existent vector of velocities (must have the right type...):
    resize!(vs, length(fs))
    for (v, f) ∈ zip(vs, fs)
        resize!(v, length(nodes(f)))  # resize velocities of a single filament
    end
    # Now we can read the velocities
    read!(io, vs, "Velocity")
end
```

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

    dataset_type = if haskey(gtop, "Lines")
        :PolyData
    elseif haskey(gtop, "Types")
        :UnstructuredGrid
    else
        error("couldn't detect the type of dataset (PolyData or UnstructuredGrid)")
    end

    gcells = if dataset_type == :PolyData
        gcells = HDF5.open_group(gtop, "Lines")  # cell info is in /VTKHDF/Lines
    elseif dataset_type == :UnstructuredGrid
        gcells = gtop  # cell info is in /VTKHDF
    end

    periods = let Ls = Vector{T}(undef, 3)
        if haskey(gtop, "Periods")
            # Files generated by VortexPasta ≥ 0.18.
            let dset = HDF5.open_dataset(gtop, "Periods")
                @assert length(dset) == length(Ls)
                read_dataset!(Ls, dset)
                close(dset)
            end
        else
            # Old versions (before VortexPasta 0.18)
            fill!(Ls, 0)  # 0 means no periodicity
        end
        ntuple(i -> Ls[i], Val(3))
    end

    with_periods = any(!iszero, periods)

    dset_points = HDF5.open_dataset(gtop, "Points")
    dset_offset = HDF5.open_dataset(gcells, "Offsets")

    Ñ, num_points_refined = size(dset_points) :: Dims{2}
    @assert Ñ == 3
    num_cells = length(dset_offset) - 1
    num_cells > 0 || error("expected at least one filament in the file")

    if with_periods
        filament_ids = Vector{FilamentId}(undef, num_cells)
        let dset = HDF5.open_dataset(gtop, "CellData/FilamentIds")
            @assert length(dset) == num_cells
            read_dataset!(filament_ids, dset)
            close(dset)
        end
    else
        filament_ids = collect(FilamentId(1):FilamentId(num_cells)) :: Vector{FilamentId}
    end

    # Determine number of actual filaments from filament ids.
    num_filaments = filament_ids[end] - filament_ids[begin] + 1
    @assert num_cells ≥ num_filaments

    if refinement == 1
        num_points = num_points_refined
    else
        # We subtract endpoint from each cell before dividing by refinement level
        num_points, _remainder = divrem(num_points_refined - num_filaments, refinement)
        @assert _remainder == 0      # the division is exact
        num_points += num_filaments  # reinclude the endpoints in final count
    end

    buf = Bumper.default_buffer()
    fs = @no_escape buf let
        # Read points as an AbstractArray{T} even if it was written as an AbstractArray{S} with S ≠ T
        # (for instance, if data was written as S = Float64 and we want to read it as T = Float32).
        # Similar thing for offsets. This ensures type stability...
        points = @alloc(Vec3{T}, num_points)
        cell_offsets = @alloc(Int, num_cells + 1)
        filament_offsets = @alloc(Int, num_filaments + 1)
        periodic_offsets = @alloc(Vec3{FilamentPeriodicOffset}, num_filaments)

        num_points_connectivity =
            only(HDF5.read_dataset(gcells, "NumberOfConnectivityIds")::Vector{Int})
        connectivity = @alloc(Int, num_points_connectivity)
        let dset = HDF5.open_dataset(gcells, "Connectivity")
            read_dataset!(connectivity, dset)
        end

        if with_periods
            let dset = HDF5.open_dataset(gtop, "FilamentPeriodicOffsets")
                memspace = HDF5.dataspace(dset)
                memtype = HDF5.datatype(FilamentPeriodicOffset)
                read_dataset!(periodic_offsets, dset; memspace, memtype)
                close(memspace)
                close(memtype)
                close(dset)
            end
        else
            fill!(periodic_offsets, zero(eltype(periodic_offsets)))
        end

        read_dataset!(cell_offsets, dset_offset)

        # Modify cell_offsets if there are gaps in the connectivity vector.
        # This corresponds to locations where the endpoint of a filament was not included in
        # the connectivity vector (to avoid ugly jumps in visualisations).
        # This must be done before compute_filament_offsets!.
        # We can ignore this if periodic wrapping was not performed in the VTKHDF file.
        if with_periods
            ngaps = 0
            for i ∈ eachindex(cell_offsets)[2:end-1]
                off = cell_offsets[i]
                a = connectivity[off]      # last index of cell i - 1
                b = connectivity[off + 1]  # first index of cell i
                Δ = b - a - 1
                @assert Δ ∈ (0, 1)
                if Δ == 1  # there's a gap! (endpoint was excluded from connectivity)
                    ngaps += 1
                end
                cell_offsets[i] += ngaps
            end
            # Check for gap after the last cell.
            # If there is no gap, then (cell_offsets[end] + ngaps) == num_points_refined.
            let off = cell_offsets[end]
                Δ = num_points_refined - (off + ngaps)
                @assert Δ ∈ (0, 1)
                if Δ == 1  # there's a gap!
                    ngaps += 1
                end
                cell_offsets[end] += ngaps
            end
        end

        compute_filament_offsets!(filament_offsets, cell_offsets, filament_ids, refinement)
        @assert last(filament_offsets) == num_points

        # Read locations in file, skipping interpolated points if refinement > 1.
        if refinement == 1
            memspace = HDF5.dataspace((Ñ, num_points))
            memtype = HDF5.datatype(T)
            read_dataset!(points, dset_points; memspace, memtype)
            close(memspace)
            close(memtype)
        else
            dspace = HDF5.dataspace(dset_points)
            n = 0
            for i ∈ 1:num_filaments
                a = filament_offsets[i] + 1
                b = filament_offsets[i + 1]
                out = @view points[a:b]
                Np = length(out)
                inds_file = range(n + 1; length = Np, step = refinement)
                HDF5.select_hyperslab!(dspace, (1:Ñ, inds_file))
                memspace = HDF5.dataspace((Ñ, Np))
                memtype = HDF5.datatype(T)
                read_dataset!(out, dset_points, dspace; memspace, memtype)
                close(memtype)
                close(memspace)
                n = last(inds_file)
            end
            close(dspace)
        end

        # Reconstruct all filaments
        map(eachindex(periodic_offsets)) do i
            periodic_offset = periodic_offsets[i]
            _load_filament(points, filament_offsets, method, i, periodic_offset, periods)
        end
    end  # @no_escape

    close(dset_points)
    close(dset_offset)

    # Load and use parametrisation values (knots `ts`) if they are included in the file.
    if haskey(gtop, "PointData/Parametrisation")
        # Read parametrisation values ts[begin:(end + 1)] for each filament.
        ts_all = map(fs) do f
            ts = knots(f)
            @view ts[begin:(end + 1)]
        end
        let dset = HDF5.open_dataset(gtop, "PointData/Parametrisation")
            _read_data_on_filaments!(T, dset, ts_all, refinement; endpoint = true)
            close(dset)
        end
        # Compute interpolation and derivative coefficients, but preserving the input knots.
        foreach(fs) do f
            ts = knots(f)
            Tper = ts[end + 1] - ts[begin]
            pad_periodic!(ts, Tper)  # apply padding to remaining parametrisation values
            if Filaments.check_nodes(Bool, f)
                update_coefficients!(f; knots = ts)
            end
        end
    else
        # Compute interpolation and derivative coefficients, also recomputing knots.
        foreach(fs) do f
            if Filaments.check_nodes(Bool, f)
                update_coefficients!(f)
            end
        end
    end

    number_of_filaments_in_file = length(fs)

    # Remove filaments that don't have enough nodes to be represented by the discretisation
    # method. This can happen e.g. if the files were written with filaments using
    # CubicSplineMethod (requires ≥ 3 nodes), and loaded using QuinticSplineMethod (requires ≥ 5 nodes).
    for i in reverse(eachindex(fs))
        if !Filaments.check_nodes(Bool, fs[i])
            popat!(fs, i)
        end
    end

    if length(fs) < number_of_filaments_in_file
        @warn(
            "Loaded number of filaments is smaller than number of filaments in file. This can happen if the discretisation method changed.",
            HDF5.filename(io), number_of_filaments_in_file, length(fs), method,
        )
    end

    VTKHDFFile(gtop, fs, refinement)
end

function compute_filament_offsets!(filament_offsets, cell_offsets, filament_ids, refinement)
    Base.require_one_based_indexing(cell_offsets)
    Base.require_one_based_indexing(filament_ids)
    num_filaments = filament_ids[end] - filament_ids[begin] + 1
    num_cells = length(cell_offsets) - 1
    @assert length(filament_ids) == num_cells ≥ num_filaments
    @assert length(filament_offsets) == num_filaments + 1
    filament_offsets[begin] = 0
    icell = firstindex(cell_offsets)
    off_i = cell_offsets[icell]  # offset of previous cell
    for i ∈ eachindex(filament_offsets)[1:end-1]
        # Count the number of *refined* points in the filament i.
        id = filament_ids[icell]  # id of the current filament
        icell += 1
        off = cell_offsets[icell]
        Np_refined = off - off_i
        off_i = off
        # If the filament was broken into multiple cells, we iterate through the other cells
        # as well (which have the same filament id).
        while icell ≤ lastindex(filament_ids) && filament_ids[icell] == id
            icell += 1
            off = cell_offsets[icell]
            Np_refined += off - off_i
            off_i = off
        end
        Np, _remainder = divrem(Np_refined - 1, refinement)
        @assert _remainder == 0
        Np += 1  # reinclude the endpoint in final count
        filament_offsets[i + 1] = filament_offsets[i] + Np
    end
    @assert icell == lastindex(cell_offsets)  # we reached the last cell
    filament_offsets
end

periodic_offset_to_apply(x⃗, x⃗_prev, Ls::Tuple, Lhs::Tuple) =
    oftype(x⃗, map(periodic_offset_to_apply, x⃗, x⃗_prev, Ls, Lhs))

function periodic_offset_to_apply(x, x_prev, L::Real, Lh::Real)
    iszero(L) && return zero(L)  # no offset if the period is "zero" (meaning no periodicity here)
    δx = x - x_prev  # we assume |δx| < L = 2 * Lh
    if δx ≥ Lh
        -L
    elseif δx ≤ -Lh
        +L
    else
        zero(L)
    end
end

# Load filament from the list of points using the given offsets.
# TODO reuse VTK cell information to detect jumps? (should be slightly faster)
function _load_filament(
        Xs::AbstractVector{Vec3{T}},
        filament_offsets::AbstractVector{<:Integer},
        method,
        i::Integer,                        # filament index
        periodic_offset::Vec3{<:Integer},  # periodic offset of first point (integer value)
        periods::NTuple{3},                # real values (0 for no periodicity)
    ) where {T}
    a = filament_offsets[i]      # Xs[a + 1] is the first point of this filament
    b = filament_offsets[i + 1]  # Xs[b] is the filament endpoint
    num_nodes = b - a - 1        # not counting the endpoint
    f = Filaments.init(ClosedFilament{T}, num_nodes, method)
    Xoffset = zero(eltype(Xs))
    current_offset = periodic_offset .* Vec3(periods)  # offset in distance units
    s⃗_prev = Xs[a + 1]
    Lhs = map(L -> L / 2, periods)
    @inbounds for j ∈ eachindex(f)
        s⃗ = Xs[a + j]  # this point is always inside the main periodic box
        p⃗ = periodic_offset_to_apply(s⃗, s⃗_prev, periods, Lhs)  # relative offset wrt previous point
        current_offset = current_offset + p⃗
        s⃗_prev = s⃗
        f[j] = s⃗ + current_offset
    end
    # Determine end-to-end offset from the endpoint.
    s⃗ = Xs[b]
    current_offset = current_offset + periodic_offset_to_apply(s⃗, s⃗_prev, periods, Lhs)
    x⃗_end = s⃗ + current_offset
    Xoffset = x⃗_end - f[begin]
    Filaments.change_offset(f, Xoffset)
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
    (; gtop, refinement,) = reader
    gdata = HDF5.open_group(gtop, "PointData")
    dset = HDF5.open_dataset(gdata, name)
    V = eltype(eltype(vs))  # typically <:AbstractFloat or SVector{3, <:AbstractFloat}
    T = eltype(V)
    if T === V  # scalar data (note: eltype(Float64) = Float64)
        ldims = ()
    else
        ldims = size(V)
    end
    _read_data_on_filaments!(T, dset, vs, refinement; ldims)
    close(dset)
    close(gdata)
    vs
end

# Read data on filament nodes.
function _read_data_on_filaments!(
        ::Type{T}, dset, vs::AbstractVector{<:AbstractVector}, refinement;
        ldims::Dims = (),
        endpoint = false,  # if `true`, it requires eltype(vs) <: PaddedVector
    ) where {T <: Number}
    V = eltype(eltype(vs))
    @assert T === eltype(V)
    N = 1 + length(ldims)
    ldims_dset..., num_points = size(dset) :: Dims{N}
    @assert ldims == ldims_dset
    n = 0
    dspace = HDF5.dataspace(dset)
    lranges = map(N -> 1:N, ldims)
    for v ∈ vs
        Np = length(v)  # number of filament nodes (may or may not include the endpoint)
        memspace = HDF5.dataspace((ldims_dset..., Np))
        memtype = HDF5.datatype(T)
        if refinement == 1
            HDF5.select_hyperslab!(dspace, (lranges..., (n + 1):(n + Np),))  # read part of the dataset
            read_dataset!(v, dset, dspace; memspace, memtype)
            n += Np + (!endpoint)  # the + (!endpoint) is because the endpoint is always included in the file (but maybe ignored here)
        else
            inds_file = range(n + 1; length = Np, step = refinement)
            HDF5.select_hyperslab!(dspace, (lranges..., inds_file,))
            read_dataset!(v, dset, dspace; memspace, memtype)
            n = last(inds_file) + refinement * (!endpoint)
        end
        maybe_pad_periodic!(v)
    end
    @assert n == num_points
    nothing
end

# Periodically pad output vector if it's a PaddedVector.
maybe_pad_periodic!(v::PaddedVector) = pad_periodic!(v)
maybe_pad_periodic!(v::AbstractVector) = v

function read_dataset!(
        out::Union{AbstractArray, Ptr},
        dset::HDF5.Dataset,
        dspace = HDF5.dataspace(dset);  # this can also be a hyperslab
        memspace = HDF5.dataspace(out),
        memtype = HDF5.datatype(out),
    )
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
