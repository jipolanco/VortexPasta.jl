using VortexPasta.BasicTypes: VectorOfVectors
using JET: @test_opt
using Test

function broadcast_factorised!(ws, us, vs)
    fill!(ws, 0)
    @. ws = 3 * (us + vs) + 2 * us
    ws
end

function broadcast_sum(us, vs)
    us .+ vs
end

@testset "VectorOfVectors" begin
    ##
    data = [rand(n) for n ∈ 1:4]
    us = @inferred VectorOfVectors(data)
    @test data === @inferred parent(us)

    @test length(us) == 4
    @test eachindex(us) === Base.OneTo(4)
    @test eachindex(us, us) === Base.OneTo(4)
    @test eltype(us) === typeof(data[1])
    @test IndexStyle(us) === IndexLinear()
    @test Base.IteratorSize(us) === Base.HasLength()
    @test Base.IteratorEltype(us) === Base.HasEltype()
    @test isdefined(similar(us), 1)

    @testset "convert" begin
        # Converting a VectorOfVectors onto a VectorOfVectors simply returns the same object.
        @test us === @inferred convert(VectorOfVectors, us)
        # Converting another kind of vector to VectorOfVectors preserves the data (doesn't
        # make a copy).
        @test us === @inferred convert(VectorOfVectors, data)
        @test data === parent(convert(VectorOfVectors, data))
    end

    @testset "push! | pop! | popat!" begin
        ws = copy(us)
        wlast = rand(5)
        push!(ws, wlast)
        @test length(ws) == length(us) + 1
        @test @views ws[1:end - 1] == us
        @test ws[end] === wlast

        w = pop!(ws)
        @test w === wlast
        @test length(ws) == length(us)
        @test ws == us

        w = popat!(ws, 2)
        @test w == us[2]
        @test length(ws) == length(us) - 1
        @test ws[1] == us[1]
        @test @views ws[2:end] == us[3:end]
    end

    @testset "copy" begin
        vs = copy(us)
        @test typeof(vs) === typeof(us)
        @test vs == us
        @test vs[1] !== us[1]  # data is not aliased, vectors were recursively copied
    end

    # `map` should return a VectorOfVectors if the given function maps vector → vector.
    # If instead it maps vector → scalar, it should return a simple vector of scalars.
    @testset "map" begin
        vs = @inferred map(copy, us)  # `copy` maps vector → vector
        vs_alt = VectorOfVectors(map(copy, parent(us)))
        @test vs == vs_alt
        @test typeof(vs) == typeof(us)  # vs is a VectorOfVectors
        @test typeof(vs) == typeof(vs_alt)
        ws = @inferred map(sum, us)  # `sum` maps vector → scalar
        ws_alt = map(sum, parent(us))
        @test ws == ws_alt
        @test typeof(ws) === typeof(ws_alt)
    end

    @testset "Broadcasting" begin
        @test_opt broadcast_sum(us, us)
        vs = @inferred us .+ 2 .* us
        @test typeof(vs) === typeof(us)
        @test vs ≈ 3 * us
        ws = similar(vs)
        @test typeof(ws) === typeof(us)
        @. ws = vs + us
        @test ws == vs + us

        # Check that there are no inference issues and no unwanted allocations.
        @test_opt broadcast_factorised!(ws, us, vs)
        broadcast_factorised!(ws, us, vs)  # run once just to be sure that everything is compiled
        @test 0 == @allocated broadcast_factorised!(ws, us, vs)
    end
    ##
end
