# Set list of points to be written to VTKHDF file.
function set_points!(
        Xs_all::AbstractVector,  # this is generally a vector of Vec3
        offsets::AbstractMatrix{<:Integer},  # (3, num_filaments)
        include_endpoint_in_connectivity::AbstractVector{Bool},
        fs, refinement, Ls::Tuple;
    )
    Base.require_one_based_indexing(offsets)
    Base.require_one_based_indexing(include_endpoint_in_connectivity)
    @assert size(offsets) == (3, length(fs))
    @assert length(include_endpoint_in_connectivity) == length(fs)
    n = firstindex(Xs_all)
    filament_ids = FilamentId[]  # this will have length Nc (= number of VTK cells)
    cell_offsets = [0]           # this will have length Nc + 1
    sizehint!(filament_ids, length(fs))
    sizehint!(cell_offsets, length(fs) + 1)
    for (nf, f) ∈ enumerate(fs)
        nprev = n
        Xs = nodes(f)
        ilast = lastindex(Xs)
        ifil = (firstindex(Xs), one(refinement))
        # Create one or more cells to describe the filament.
        while ifil[1] ≤ ilast
            n_new, ifil = add_filament_cell!(Xs_all, n, f, refinement, Ls, ifil)
            n_written = n_new - n
            n = n_new
            push!(filament_ids, nf)  # filament id associated to created cell
            push!(cell_offsets, cell_offsets[end] + n_written)
        end
        # Ignore the endpoint if there is a jump between the last 2 points of the filament
        # (where the last point is the endpoint f[end + 1]).
        include_endpoint = Bool(ifil[2])
        include_endpoint_in_connectivity[nf] = include_endpoint
        if !include_endpoint
            cell_offsets[end] -= 1
        end
        @assert n == nprev + length(Xs) * refinement + 1
        # Look at the periodic offset of the first point of the filament.
        # We actually look at the *last* written point (Xs_all[n - 1]) for memory locality
        # (performance) reasons. Note that this point is always equal to Xs_all[nprev] (the
        # first written point of the filament).
        let δx⃗ = Xs[begin] - Xs_all[n - 1]
            I = eltype(offsets)
            p⃗ = map(δx⃗, Ls) do δx, L
                L === nothing && return zero(I)
                p = round(I, δx / L)
                # If p == 0, one can get a δx which is not exactly zero (but of the order of 1e-16).
                # @assert (iszero(p) && abs(δx) < 10 * eps(L)) || (p * L ≈ δx)
                p
            end
            for (i, p) ∈ pairs(p⃗)
                offsets[i, nf] = p
            end
        end
    end
    @assert n == lastindex(Xs_all) + 1
    (; filament_ids, cell_offsets,)
end

# Create the specification of a single VTK cell. A cell stops either at the end of a
# filament, or when the filament crosses the domain boundaries (whichever happens first).
function add_filament_cell!(
        Xs_all, n::Int, f::ClosedFilament,
        refinement, Ls,
        ifil = (firstindex(nodes(f)), 1),
    )
    Xs = nodes(f)
    i, j = ifil
    @assert n ∈ eachindex(Xs_all)  # Xs_all[n] is the next point to be set
    @assert i ∈ eachindex(Xs)
    @assert j ∈ 1:refinement
    Lhs = map(L -> L === nothing ? nothing : L / 2, Ls)  # half periods
    ξs = range(0, 1; length = refinement + 1)[1:refinement]
    s⃗_prev = let x⃗ = j == 1 ? f[i] : f(i, ξs[j])
        Filaments.to_main_periodic_cell(x⃗, Ls)
    end
    while i ≤ lastindex(Xs)
        while j ≤ refinement
            x⃗ = j == 1 ? f[i] : f(i, ξs[j])
            s⃗ = Filaments.to_main_periodic_cell(x⃗, Ls)
            if Filaments.is_jump(s⃗, s⃗_prev, Lhs)
                # The filament (likely) crossed a domain boundary.
                # We stop here and recompute this point when creating the next cell.
                return n, (i, j)
            end
            Xs_all[n] = s⃗
            s⃗_prev = s⃗
            n += 1
            j += 1
        end
        i += 1
        j = 1
    end
    @assert i == lastindex(Xs) + 1
    @assert j == 1
    # Close the curve
    s⃗ = Filaments.to_main_periodic_cell(Xs[i], Ls)
    Xs_all[n] = s⃗
    n += 1
    # Return j = 0 if the endpoint is "far" from the previous point after periodic wrapping.
    # In that case, we won't include the endpoint in the cell connectivity vector.
    # This is just done to avoid ugly jumps in visualisations.
    j = Filaments.is_jump(s⃗, s⃗_prev, Lhs) ? 0 : 1
    n, (i, j)  # here i == lastindex(Xs) + 1 and j ∈ {0, 1}
end

# Write point locations to VTKHDF dataset.
# We use the low-level HDF5 API since it allows us to efficiently reinterpret Xs_all as a
# (3, Np) matrix of real values (as it works directly with pointers).
function points_to_vtkhdf_dataset(gtop, Xs_all::AbstractVector)
    V = eltype(Xs_all)  # usually SVector{3, T}
    N = length(V)       # usually == 3 (works when V <: StaticArray)
    @assert N == 3
    T = eltype(V)
    @assert T <: AbstractFloat
    num_points = length(Xs_all)
    dtype = HDF5.datatype(T)
    dspace = HDF5.dataspace((N, num_points))
    dset = HDF5.create_dataset(gtop, "Points", dtype, dspace)
    memtype = dtype
    memspace = dspace
    HDF5.API.h5d_write(dset.id, memtype.id, memspace.id, dspace.id, dset.xfer, Xs_all)
    close(dset)
    close(dspace)
    close(dtype)
    nothing
end
