using VortexPasta.PaddedArrays
using VortexPasta.PaddedArrays: FromRight, FromLeft, FromCentre, npad
using Test

@testset "PaddedArrays" begin
    @testset "Basic functionality" begin
        data = collect(1:20)
        v = @inferred PaddedArray{2}(data)
        M = npad(v)
        @test M == 2
        N = length(v)
        @test N == length(data) - 2M
        @test eachindex(v) == 1:N
        for i ∈ 1:M
            # Check that we can access ghost values
            @test checkbounds(Bool, v, 1 - i)
            @test checkbounds(Bool, v, N + i)
        end
        @test !checkbounds(Bool, v, -M)
        @test !checkbounds(Bool, v, N + M + 1)
        @test v == @view data[(M + 1):(M + N)]
    end

    @testset "Periodic padding" begin
        data = collect(1:20)
        v = @inferred PaddedArray{2}(data)
        M = npad(v)
        @test M == 2
        vc = collect(v) :: Vector  # create copy of central data

        pad_periodic!(v)
        @test vc == v  # central values didn't change
        # Check that ghost values have been correctly modified
        @test v[begin - 2] == vc[end - 1]
        @test v[begin - 1] == vc[end]
        @test v[end + 1] == vc[begin]
        @test v[end + 2] == vc[begin + 1]

        L = 42
        pad_periodic!(v, L)
        @test v[begin - 2] == vc[end - 1] - L
        @test v[begin - 1] == vc[end] - L
        @test v[end + 1] == vc[begin] + L
        @test v[end + 2] == vc[begin + 1] + L
    end

    @testset "Periodic padding: FromLeft" begin
        data = collect(1:20)
        v = @inferred PaddedArray{2}(copy(data))
        M = npad(v)
        @test M == 2
        vc = collect(v) :: Vector  # create copy of central data
        pad_periodic!(FromLeft(), v)
        @test vc != v  # some central values changed
        # Check that ghost values have been correctly modified
        @test v[begin - 2] == v[end - 1] == data[begin]  # left data is preserved
        @test v[begin - 1] == v[end] == data[begin + 1]  # left data is preserved
        @test v[end + 1] == v[begin] == data[begin + 2]
        @test v[end + 2] == v[begin + 1] == data[begin + 3]
    end

    @testset "Periodic padding: FromRight" begin
        data = collect(1:20)
        v = @inferred PaddedArray{2}(copy(data))
        M = npad(v)
        @test M == 2
        vc = collect(v) :: Vector  # create copy of central data
        pad_periodic!(FromRight(), v)
        @test vc != v  # some central values changed
        # Check that ghost values have been correctly modified
        @test v[begin - 2] == v[end - 1] == data[end - 3]
        @test v[begin - 1] == v[end] == data[end - 2]
        @test v[end + 1] == v[begin] == data[end - 1]  # right data is preserved
        @test v[end + 2] == v[begin + 1] == data[end]  # right data is preserved
    end

    @testset "Periodic padding: edge cases" begin
        # This is a special case where the length of `v` is less than the padding `M`.
        data = collect(1:8)
        v = @inferred PaddedArray{3}(copy(data))
        N = length(v)
        M = npad(v)
        @assert N == 2 && M == 3
        L = 42
        @testset "FromCentre" begin
            u = copy(v)
            pad_periodic!(FromCentre(), u, L)
            @test u[1] == data[M + 1]  # didn't change
            @test u[2] == data[M + 2]  # didn't change
            @test u[0] == u[2] - L
            @test u[-1] == u[1] - L
            @test u[-2] == u[2] - 2L
            @test u[N + 1] == u[1] + L
            @test u[N + 2] == u[2] + L
            @test u[N + 3] == u[1] + 2L
        end
        @testset "FromLeft" begin
            u = copy(v)
            pad_periodic!(FromLeft(), u, L)
            @test u[0] == data[M - 0]   # didn't change
            @test u[-1] == data[M - 1]  # didn't change
            @test u[-2] == u[0] - L
            for i ∈ 1:(N + M)
                @test u[i] == u[i - 2] + L
            end
        end
        @testset "FromRight" begin
            u = copy(v)
            pad_periodic!(FromRight(), u, L)
            @test u[N + 1] == data[M + N + 1]  # didn't change
            @test u[N + 2] == data[M + N + 2]  # didn't change
            @test u[N + 3] == u[N + 1] + L
            for i ∈ (1 - M):N
                @test u[i] == u[i + 2] - L
            end
        end
    end

    @testset "Views" begin
        data = collect(1:20)
        v = @inferred PaddedArray{2}(data)
        let w = @inferred view(v, :)
            @test w isa SubArray
            @test parent(w) == data  # this is a view to `data`, not to `v`
            @test length(w) == length(v)  # ghost cells are not included
            @test w == v
        end
        let w = @inferred view(v, 2:3)
            @test w isa SubArray
            @test parent(w) == data  # this is a view to `data`, not to `v`
            @test length(w) == 2
            @test w == v[2:3]
        end
        let w = @inferred view(v, 3)
            @test w isa SubArray
            @test parent(w) == data  # this is a view to `data`, not to `v`
            @test ndims(w) == 0
            @test w[] == v[3]
        end
    end

    @testset "Equality" begin
        data = collect(1:20)
        v = @inferred PaddedArray{2}(data)
        w = similar(v)
        fill!(w, 0)
        # 1. Copy only non-ghost cells.
        @views w[:] .= v[:]  # this doesn't include ghost cells
        @test w != v  # arrays are different because ghost cells were not copied!
        @test !isapprox(w, v)
        # 2. Copy including ghost cells.
        copyto!(w, v)  # this includes ghost cells
        @test w == v   # arrays are now equal
        @test isapprox(w, v)
    end
end
