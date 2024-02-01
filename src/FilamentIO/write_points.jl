function write_points(points, fs, refinement)
    if refinement == 1
        write_points_unrefined(points, fs)
    else
        write_points_refined(points, fs, refinement)
    end :: Int
end

function write_points_unrefined(points, fs)
    V = eltype(eltype(fs))  # usually V == SVector{3, T}
    N = length(V)  # usually == 3 (works when V <: StaticArray)
    n = 0
    for f ∈ fs
        Xs = nodes(f)
        Np = length(Xs) + 1  # include the endpoint
        write_padded_data_to_hdf5(points, Xs, n, (N,), Np)  # write data to existing HDF5 dataset
        n += Np
    end
    n
end

function write_points_refined(points, fs, refinement)
    V = eltype(eltype(fs))  # usually V == SVector{3, T}
    N = length(V)  # usually == 3 (works when V <: StaticArray)
    n = 0
    buf = Bumper.default_buffer()
    for f ∈ fs
        @no_escape buf begin
            Np = refinement * length(f) + 1  # number of output points (the +1 is to close the loop)
            Xs = @alloc(V, Np)
            refine_filament_coordinates!(Xs, f, refinement)
            @assert Xs[end] == f[end + 1]  # already includes the endpoint
            write_padded_data_to_hdf5(points, Xs, n, (N,), Np)  # write data to existing HDF5 dataset
            n += Np
        end
    end
    n
end

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

function refine_filament_coordinates!(Xs::AbstractVector, f::ClosedFilament, refinement::Int)
    refinement ≥ 1 || error("refinement must be ≥ 1")
    n = firstindex(Xs, 2) - 1
    for j ∈ eachindex(f)
        Xs[n += 1] = f[j]
        for m ∈ 2:refinement
            ζ = on_segment_location(m, refinement)  # in ]0, 1[
            Xs[n += 1] = f(j, ζ)
        end
    end
    @assert n == refinement * length(f)
    Xs[n += 1] = f[end + 1]  # close the loop
    Xs
end
