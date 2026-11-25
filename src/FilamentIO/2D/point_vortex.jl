using ..BiotSavart: BiotSavart, PointVortexState

# Note: this file reuses definitions in ../VTKHDF/vtkhdf.jl

struct PointVortexVTKHDFFile{
        Group <: HDF5.Group,
        Vortices <: PointVortexState,
    }
    gtop::Group  # "/VTKHDF" group
    vortices::Vortices
end

root(io::PointVortexVTKHDFFile) = HDF5.root(io.gtop)::HDF5.Group  # returns the "/" group

function write_vtkhdf(
        func::Func,
        filename::AbstractString,
        vortices::PointVortexState;
        periods = nothing,
    ) where {Func <: Function}
    if periods isa Tuple
        periods_ = periods
    else
        @assert periods isa Union{Nothing, Real}
        periods_ = (periods, periods)  # same period in all directions (works with `nothing` too)
    end
    HDF5.h5open(filename, "w") do io
        writer = init_vtkhdf(io, vortices; periods = periods_)
        func(writer)
    end
    nothing
end

function init_vtkhdf(
        io::HDF5.File,
        vortices::PointVortexState;
        periods::Tuple,
    )
    gtop = HDF5.create_group(io, "VTKHDF")
    HDF5.attrs(gtop)["Version"] = [2, 1]
    HDF5.attrs(gtop)["VortexPasta minimum version"] = [0, 36, 2]

    (; points, signs) = vortices
    V = eltype(points)  # usually V == SVector{2, T}
    T = eltype(V)

    gtop["Periods"] = collect(
        map(periods) do L
            L === nothing ? zero(T) : T(L)
        end
    )

    # This attribute *must* be written as ASCII instead of UTF8, or ParaView will fail to
    # open the file.
    dataset_type = :PolyData
    let s = string(dataset_type)
        dtype = datatype_ascii(s)
        dspace = HDF5.dataspace(s)
        attr = HDF5.create_attribute(gtop, "Type", dtype, dspace)
        HDF5.write_attribute(attr, dtype, s)
    end

    # Add empty categories (note: we have no cells at all)
    if dataset_type == :PolyData
        for type ∈ ("Lines", "Polygons", "Strips", "Vertices")
            gempty = HDF5.create_group(gtop, type)
            gempty["NumberOfCells"] = [0]
            gempty["NumberOfConnectivityIds"] = [0]
            gempty["Connectivity"] = Int[]
            gempty["Offsets"] = [0]
            close(gempty)
        end
    end

    num_points = length(points)
    gtop["NumberOfPoints"] = [num_points]

    g_pointdata = HDF5.create_group(gtop, "PointData")
    g_celldata = HDF5.create_group(gtop, "CellData")

    vector_data_to_vtkhdf_dataset(gtop, points)  # write "/VTKHDF/Points" dataset
    g_pointdata["signs"] = one(Int8) .* signs    # this converts signs (±1) to Int8

    close(g_celldata)
    close(g_pointdata)

    PointVortexVTKHDFFile(gtop, vortices)
end

function Base.setindex!(
        writer::PointVortexVTKHDFFile,
        vs::AbstractVector,
        name::AbstractString,
    )
    write_point_data(writer, vs, name)
end

function write_point_data(writer::PointVortexVTKHDFFile, vs, name)
    (; gtop,) = writer
    gdata = HDF5.open_group(gtop, "PointData")
    T = eltype(vs)  # usually <:Real or Vec3{<:Real}
    attrname = vtk_attribute_type(T)  # "Scalars", "Vectors", ...
    names = if haskey(HDF5.attrs(gdata), attrname)
        names_prev = HDF5.attrs(gdata)[attrname]::Vector{String}
        vcat(names_prev, name)::Vector{String}
    else
        String[name]
    end
    HDF5.attrs(gdata)[attrname] = names  # overwrite attribute if it already existed
    write_point_data(gdata, name, vs)
    close(gdata)
    nothing
end
