using VortexPasta.PaddedArrays: PaddedVector, pad_periodic!
using VortexPasta.Filaments
using VortexPasta.FilamentIO
using VortexPasta.BasicTypes: VectorOfVectors
using VortexPasta.BiotSavart
using LinearAlgebra: ⋅
using HDF5: HDF5
using Test

function init_ring_filament(; R, z, sign)
    tlims = (0.0, 2.0)
    S(t) = Vec3(π + R * cospi(sign * t), π + R * sinpi(sign * t), z)
    (; R, z, sign, tlims, S,)
end

function tangent_streamfunction!(ψt_all, fs, ψs_all::VectorOfVectors)
    for (f, ψs, ψt) ∈ zip(fs, ψs_all, ψt_all)
        for i ∈ eachindex(f, ψs, ψt)
            ψt[i] = ψs[i] ⋅ f[i, UnitTangent()]
        end
    end
    ψt_all
end

function test_hdf5()
    # Copied from the ring collision test.
    R = π / 3  # ring radius
    L = π / 8  # ring distance
    rings = [
        init_ring_filament(; R, z = π - L / 2, sign = +1),
        init_ring_filament(; R, z = π + L / 2, sign = -1),
    ]
    fs = map(rings) do ring
        (; tlims, S,) = ring
        N = 64
        ζs = range(tlims...; length = N + 1)[1:N]
        Filaments.init(ClosedFilament, S.(ζs), CubicSplineMethod())
    end

    Γ = 2.0
    a = 1e-6
    Δ = 1/4
    Ls = (2π, 2π, 2π)
    Ns = (64, 64, 64)
    kmax = (Ns[1] ÷ 2) * 2π / Ls[1]
    α = kmax / 6
    rcut = 4 / α

    params = ParamsBiotSavart(;
        Γ, α, a, Δ, rcut, Ls, Ns,
    )

    cache = @inferred BiotSavart.init_cache(params, fs)

    # Use VectorOfVectors to make things more interesting (this is not necessary, but can be
    # convenient)
    vs = VectorOfVectors(map(similar ∘ nodes, fs))
    ψs = similar(vs)
    fields = (velocity = vs, streamfunction = ψs)
    @inferred BiotSavart.compute_on_nodes!(fields, cache, fs)

    ψt = similar(ψs, Float64) :: AbstractVector{<:PaddedVector}
    tangent_streamfunction!(ψt, fs, ψs)

    foreach(pad_periodic!, vs)
    foreach(pad_periodic!, ψs)
    foreach(pad_periodic!, ψt)

    # Try writing data which is *not* backed by a PaddedVector.
    # When refinement is enabled, values on the last segment should be obtained by
    # interpolation using the last and first data points, since there's no padding.
    ψt_vec = VectorOfVectors(map(ψ -> Vector{Float64}(undef, length(ψ)), ψs)) :: AbstractVector{<:Vector}
    tangent_streamfunction!(ψt_vec, fs, ψs)
    @test ψt == ψt_vec

    # Curvature along filaments
    curvature_vector = VectorOfVectors(map(fs) do f
        map(i -> f[i, CurvatureVector()], eachindex(f))
    end)
    curvature_scalar = VectorOfVectors(map(fs) do f
        map(i -> f[i, CurvatureScalar()], eachindex(f))
    end)

    time = 0.3
    info_str = ["one", "two"]
    dataset_types = (:PolyData, :UnstructuredGrid)

    function check_fields(io)
        vs_read = @inferred read(io, "velocity", PointData(), Vec3{Float64})
        ψs_read = @inferred read(io, "streamfunction", PointData(), Vec3{Float64})
        ψt_read = @inferred read(io, "streamfunction_t", PointData(), Float64)

        @test vs == vs_read
        @test ψs == ψs_read
        @test ψt == ψt_read

        # Check that arrays are correctly padded
        for (us, vs) ∈ (vs => vs_read, ψs => ψs_read, ψt => ψt_read)
            for (u, v) ∈ zip(us, vs)
                @test v isa PaddedVector
                @test parent(v) == parent(u)  # this also compares "ghost" entries
            end
        end

        # Test reading onto VectorOfVectors
        ψs_alt = similar(ψs)
        @assert ψs_alt != ψs
        @assert ψs_alt isa VectorOfVectors
        read!(io, ψs_alt, "streamfunction")
        @test ψs_alt == ψs
        for (u, v) ∈ zip(ψs, ψs_alt)
            @test v isa PaddedVector
            @test parent(v) == parent(u)  # this also compares "ghost" entries
        end

        # Test reading onto non-PaddedVectors
        ψt_vec_read = similar(ψt_vec)
        @assert ψt_vec_read isa VectorOfVectors
        @assert eltype(ψt_vec_read) <: Vector
        read!(io, ψt_vec_read, "streamfunction_t_vec")
        @test ψt_vec_read == ψt_vec

        # Test reading geometric quantities
        let us = similar(curvature_vector)
            read!(io, us, "curvature_vector")
            @test us == curvature_vector
        end
        let us = similar(curvature_scalar)
            read!(io, us, "curvature_scalar")
            @test us == curvature_scalar
        end

        # Test reading field data
        time_read = @inferred read(io, "time", FieldData(), Float64)  # this is a vector!
        @test time == only(time_read)

        info_read = read(io, "info", FieldData(), String)
        @test info_str == info_read
    end

    @testset "→ refinement = $refinement" for refinement ∈ (1, 3)
        @testset "Dataset type: $dataset_type" for dataset_type ∈ dataset_types
            fname = "ring_collision_ref$(refinement)_$dataset_type.vtkhdf"

            # Write results
            FilamentIO.write_vtkhdf(fname, fs; refinement, dataset_type, parametrisation = false) do io
                io["velocity"] = vs
                io["streamfunction"] = ψs
                io["streamfunction_t"] = ψt
                io["streamfunction_t_vec"] = ψt_vec
                io["curvature_vector"] = CurvatureVector()
                io["curvature_scalar"] = CurvatureScalar()
                io["time"] = time
                io["info"] = info_str
            end

            HDF5.h5open(fname, "r") do ff
                gbase = HDF5.open_group(ff, "VTKHDF")
                ref = read(gbase["RefinementLevel"]) :: Int
                @test ref == refinement
                # Check that both datasets are equal in the file.
                # This is mainly interesting when refinement > 1, so that we're checking also
                # interpolated values in-between nodes.
                local ψt_read = read(gbase["PointData/streamfunction_t"]) :: Vector{Float64}
                local ψt_vec_read = read(gbase["PointData/streamfunction_t_vec"]) :: Vector{Float64}
                @test ψt_read == ψt_vec_read
            end

            # Read results back
            fs_read = @inferred FilamentIO.read_vtkhdf(
                check_fields, fname, Float64, CubicSplineMethod(),
            )

            @test eltype(eltype(fs_read)) === Vec3{Float64}
            @test fs == fs_read

            fs_read_f32 = @inferred FilamentIO.read_vtkhdf(fname, Float32, CubicSplineMethod())
            @test eltype(eltype(fs_read_f32)) === Vec3{Float32}
            @test fs ≈ fs_read_f32

            # Same without passing a function. Also, modify parametrisation (knots) of the
            # filament, to make sure we read it back with the same parametrisation.
            gs = map(copy, fs)
            for g ∈ gs
                ts = knots(g)
                for i ∈ eachindex(ts)
                    ts[i] = sqrt(ts[i])
                end
                pad_periodic!(ts, ts[end + 1] - ts[begin])
                update_coefficients!(g; knots = ts)
            end
            FilamentIO.write_vtkhdf(fname * ".alt", gs; refinement)
            gs_read = FilamentIO.read_vtkhdf(fname * ".alt", Float64, CubicSplineMethod())
            for (f, g) ∈ zip(gs, gs_read)
                @test knots(f) == knots(g)      # parametrisation is the same!
                @test f.coefs.cs == g.coefs.cs  # interpolation coefficients are the same!
            end
            @test gs == gs_read
        end
    end

    nothing
end

@testset "FilamentIO: HDF5" begin
    test_hdf5()
end
