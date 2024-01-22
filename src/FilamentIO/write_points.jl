function write_points_unrefined(points, fs)
    V = eltype(eltype(fs))  # usually V == SVector{3, T}
    N = length(V)  # usually == 3 (works when V <: StaticArray)
    n = 0
    for f ∈ fs
        Xs = nodes(f)
        Np = length(Xs) + 1  # include the endpoint
        write_points_to_hdf5(points, Xs, n, (N,), Np)  # write data to existing HDF5 dataset
        n += Np
    end
    n
end

function write_points_refined(points, fs, refinement)
    V = eltype(eltype(fs))  # usually V == SVector{3, T}
    T = eltype(V)
    N = length(V)  # usually == 3 (works when V <: StaticArray)
    n = 0
    buf = Bumper.default_buffer()
    for f ∈ fs
        @no_escape buf begin
            Np = refinement * length(f) + 1  # number of output points (the +1 is to close the loop)
            Xs = @alloc(T, N, Np)
            refine_filament_coordinates!(Xs, f, refinement)
            @assert @view(Xs[:, end]) ≈ f[end + 1]  # already includes the endpoint
            write_points_to_hdf5(points, Xs, n, (N,), Np)  # write data to existing HDF5 dataset
            n += Np
        end
    end
    n
end

function write_points_to_hdf5(points, Xs, n, ldims, Np)
    (; dset, dspace, dtype,) = points
    memspace = HDF5.dataspace((ldims..., Np))
    memtype = dtype
    lranges = map(N -> 1:N, ldims)
    HDF5.select_hyperslab!(dspace, (lranges..., (n + 1):(n + Np)))
    HDF5.API.h5d_write(dset.id, memtype.id, memspace.id, dspace.id, dset.xfer, Xs)
    close(memspace)
    nothing
end

function refine_filament_coordinates!(Xs::AbstractMatrix, f::ClosedFilament, refinement::Int)
    Xs_nodes = nodes(f)
    refinement ≥ 1 || error("refinement must be ≥ 1")
    refinement == 1 && return Xs_nodes[begin:end + 1]
    V = eltype(f)  # usually V == SVector{3, T}
    n = firstindex(Xs, 2) - 1
    subinds = range(0, 1; length = refinement + 1)[1:refinement]
    for j ∈ eachindex(f), ζ ∈ subinds
        n += 1
        x⃗ = f(j, ζ) :: V
        for i ∈ eachindex(x⃗)
            Xs[i, n] = x⃗[i]
        end
    end
    @assert n == refinement * length(f)
    let x⃗ = f(lastindex(f), 1.0)  # close the loop
        n += 1
        for i ∈ eachindex(x⃗)
            Xs[i, n] = x⃗[i]
        end
    end
    Xs
end
