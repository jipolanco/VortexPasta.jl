function write_parametrisation(info, fs, refinement)
    if refinement == 1
        write_parametrisation_unrefined(info, fs)
    else
        write_parametrisation_refined(info, fs, refinement)
    end :: Int
end

function write_parametrisation_unrefined(info, fs)
    n = 0
    for f ∈ fs
        ts = knots(f)
        Np = length(ts) + 1  # include the endpoint
        ldims = ()           # case of scalar data
        write_padded_data_to_hdf5(info, ts, n, ldims, Np)  # write data to existing HDF5 dataset
        n += Np
    end
    n
end

function write_parametrisation_refined(info, fs, refinement)
    n = 0
    buf = Bumper.default_buffer()
    for f ∈ fs
        @no_escape buf begin
            ts = knots(f)
            T = eltype(ts)
            Np = refinement * length(f) + 1  # number of output points (the +1 is to close the loop)
            Xs = @alloc(T, Np)
            refine_parametrisation!(Xs, ts, refinement)
            @assert Xs[end] == ts[end + 1]  # already includes the endpoint
            ldims = ()           # case of scalar data
            write_padded_data_to_hdf5(info, Xs, n, ldims, Np)  # write data to existing HDF5 dataset
            n += Np
        end
    end
    n
end

function refine_parametrisation!(Xs::AbstractVector, ts::PaddedVector, refinement::Int)
    refinement ≥ 1 || error("refinement must be ≥ 1")
    n = firstindex(Xs, 2) - 1
    for j ∈ eachindex(ts)
        ta, tb = ts[j], ts[j + 1]
        Xs[n += 1] = ta
        for m ∈ 2:refinement
            ζ = on_segment_location(m, refinement)  # in ]0, 1[
            Xs[n += 1] = (1 - ζ) * ta + ζ * tb
        end
    end
    @assert n == refinement * length(ts)
    Xs[n += 1] = ts[end + 1]  # close the loop
    Xs
end
