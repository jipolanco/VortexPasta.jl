using VortexPasta.BasicTypes: VectorOfVectors
using Test

@testset "VectorOfVectors" begin
    data = [rand(n) for n ∈ 1:4]
    us = @inferred VectorOfVectors(data)

    @test length(us) == 4
    @test eachindex(us) === Base.OneTo(4)
    @test eachindex(us, us) === Base.OneTo(4)
    @test eltype(us) === typeof(data[1])
    @test IndexStyle(us) === IndexLinear()
    @test Base.IteratorSize(us) === Base.HasLength()
    @test Base.IteratorEltype(us) === Base.HasEltype()
    @test isdefined(similar(us), 1)

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

    @testset "Broadcasting" begin
        vs = @inferred us .+ 2 .* us
        @test typeof(vs) === typeof(us)
        @test vs ≈ 3 * us
        ws = similar(vs)
        @test typeof(ws) === typeof(us)
        @. ws = vs + us
        @test ws == vs + us
    end
end
