# Write data on filament nodes -- case of scalar data (e.g. local curvature magnitude)
function write_data_on_filaments(
        ::Type{T}, writer::VTKHDFFile, gdata, vs::AbstractVector{<:AbstractVector},
        name,
    ) where {T <: Number}
    (; refinement, fs,) = writer
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
        write_data_refined(info, vs, fs, refinement; ldims)
    end
    @assert n == num_points

    map(close, info)

    nothing
end

# Write data on filament nodes -- case of vector data (e.g. velocity)
function write_data_on_filaments(
        ::Type{V}, writer::VTKHDFFile, gdata, vs::AbstractVector{<:AbstractVector},
        name,
    ) where {V <: AbstractVector}
    (; refinement, fs,) = writer
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
        write_data_refined(info, vs, fs, refinement; ldims)
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

function write_data_refined(info, vs, fs, refinement; ldims::Dims = ())
    V = eltype(eltype(vs))  # usually SVector{3, T} or just T (scalar data)
    T = eltype(V)
    @assert T <: Number
    @assert length(vs) == length(fs)  # 1 data vector per filament
    n = 0
    buf = Bumper.default_buffer()
    for (v, f) ∈ zip(vs, fs)
        @no_escape buf begin
            Np = refinement * length(v) + 1  # number of output points (the +1 is to close the loop)
            Xs = @alloc(V, Np)
            refine_data_on_filament!(Xs, v, f, refinement, buf)
            write_points_to_hdf5(info, Xs, n, ldims, Np)
            n += Np
        end
    end
    n
end

function refine_data_on_filament!(
        Xs::AbstractVector{T}, vs::AbstractVector{T}, f::ClosedFilament,
        refinement::Int, buf,
    ) where {T}
    method = Filaments.discretisation_method(f)
    ts = Filaments.knots(f)
    M = Filaments.npad(method)
    Np = length(vs)
    @no_escape buf begin
        data = @alloc(T, Np + 2M)
        cs = PaddedVector{M}(data)
        nderiv = Filaments.required_derivatives(method)
        cderiv = ntuple(Val(nderiv)) do _
            local data = @alloc(T, Np + 2M)
            PaddedVector{M}(data)
        end
        coefs = Filaments.init_coefficients(method, cs, cderiv)
        Filaments.compute_coefficients!(coefs, vs, ts)
        n = firstindex(Xs)::Int - 1
        for j ∈ eachindex(vs)
            Xs[n += 1] = vs[j]  # copy value on node
            for m ∈ 2:refinement
                # Copy value interpolated in-between nodes
                ζ = (m - 1) / refinement  # in ]0, 1[
                Xs[n += 1] = Filaments.evaluate(coefs, ts, j, ζ)
            end
        end
        Xs[n += 1] = vs[begin]  # close the loop
        @assert n == lastindex(Xs)
    end
    Xs
end
