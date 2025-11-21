using VortexPasta.CellLists
using VortexPasta.BiotSavart: Vec3, CPU, PseudoGPU
using Adapt: adapt
using StaticArrays: SVector
using StableRNGs: StableRNG
using Random
using Test

@inline function deperiodise_separation(r⃗::Vec3, Ls, Ls_half)
    map(deperiodise_separation, r⃗, Vec3(Ls), Vec3(Ls_half))::typeof(r⃗)
end

@inline function deperiodise_separation(r::Real, L::Real, Lhalf::Real)
    while r > Lhalf
        r -= L
    end
    while r < -Lhalf
        r += L
    end
    r
end

function compute_interaction_naive(f::F, xp, vp, r_cut, Ls) where {F}
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    interaction = zero(eltype(vp))
    n_interactions = 0
    @inbounds for i in eachindex(xp, vp)
        x⃗ = xp[i]
        for j in eachindex(xp, vp)
            y⃗ = xp[j]
            r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            r² > r²_cut && continue
            interaction += f(vp[i], vp[j], r⃗, r²)
            n_interactions += 1
        end
    end
    (; interaction, n_interactions)
end

function construct_cell_list(r_cut, Ls; backend = CPU(), nsubdiv = Val(1))
    rs_cut = map(_ -> r_cut, Ls)
    @inferred PeriodicCellList(backend, rs_cut, Ls, nsubdiv)
end

function compute_interaction_nearby_elements(f::F, cl::PeriodicCellList, xp, vp, r_cut) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    @assert length(xp) == length(vp)
    CellLists.set_elements!(cl, xp)
    interaction = zero(eltype(vp))
    n_interactions = 0
    @inbounds for i in eachindex(xp, vp)
        x⃗ = xp[i]
        for j in CellLists.nearby_elements(cl, x⃗)
            y⃗ = xp[j]
            r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            r² > r²_cut && continue
            interaction += f(vp[i], vp[j], r⃗, r²)
            n_interactions += 1
        end
    end
    (; interaction, n_interactions)
end

# This is similar to above but with a closure.
function compute_interaction_foreach_source(f::F, cl::PeriodicCellList, xp, vp, r_cut) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    @assert length(xp) == length(vp)

    wp = similar(vp)  # this is the "influence" on a point of all neighbouring points
    n_interactions_p = similar(vp, Int)  # this is the number of points that affect a given point
    fill!(wp, 0)
    fill!(n_interactions_p, 0)

    CellLists.set_elements!(cl, xp)

    @inbounds Threads.@threads for i in eachindex(xp, vp)
        x⃗ = xp[i]
        CellLists.foreach_source(cl, x⃗) do j
            @inbounds y⃗ = xp[j]
            r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            @inbounds if r² <= r²_cut
                wp[i] += f(vp[i], vp[j], r⃗, r²)
                n_interactions_p[i] += 1
            end
        end
    end

    interaction = sum(wp)
    n_interactions = sum(n_interactions_p)

    (; interaction, n_interactions)
end

# This variant is more adapted for GPUs.
function compute_interaction_foreach_pair(
        f::F, cl::PeriodicCellList, xp, vp, r_cut;
        folded = Val(false),
        batch_size = nothing,
    ) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    @assert length(xp) == length(vp)

    wp = similar(vp)  # this is the "influence" on a point of all neighbouring points
    n_interactions_p = similar(vp, Int)  # this is the number of points that affect a given point
    fill!(wp, 0)
    fill!(n_interactions_p, 0)

    CellLists.set_elements!(cl, xp; folded)

    if batch_size === nothing
        CellLists.foreach_pair(cl, xp; folded) do x⃗, i, j
            # @inbounds x⃗ = xp[i]
            @inbounds y⃗ = xp[j]
            r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            @inbounds if r² <= r²_cut
                wp[i] += f(vp[i], vp[j], r⃗, r²)
                n_interactions_p[i] += 1
            end
        end
    else
        CellLists.foreach_pair(cl, xp; folded, batch_size) do x⃗, i, js, m
            # @inbounds x⃗ = xp[i]
            for n in 1:m
                j = js[n]
                @inbounds y⃗ = xp[j]
                r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
                r² = sum(abs2, r⃗)
                @inbounds if r² <= r²_cut
                    wp[i] += f(vp[i], vp[j], r⃗, r²)
                    n_interactions_p[i] += 1
                end
            end
        end
    end

    interaction = sum(wp)
    n_interactions = sum(n_interactions_p)

    (; interaction, n_interactions)
end

function test_cell_lists()
    T = Float64
    Ls = (2π, 3π, 4π)
    r_cut = π/4

    Np = 3000
    rng = StableRNG(42)

    # Random locations and values within the domain
    xp = rand(rng, Vec3{T}, Np)
    vp = randn(rng, T, Np)
    for i in eachindex(xp)
        xp[i] = xp[i] .* Ls
    end

    # Arbitrary interaction function
    function f_interaction(u, v, r⃗, r²)
        r_c = oftype(r², 2π / 100)  # to avoid singularity
        r = sqrt(r²)
        u * v / (r + r_c)
    end

    run_naive = compute_interaction_naive(f_interaction, xp, vp, r_cut, Ls)

    @testset "Number of subdivisions = $nsubdiv" for nsubdiv in 1:2
        cl = construct_cell_list(r_cut, Ls; nsubdiv = Val(nsubdiv))

        @testset "Using nearby_elements" begin
            @test isempty(cl.next_index)
            @test startswith("PeriodicCellList{3} with:")(repr(cl))
            run_cl = compute_interaction_nearby_elements(f_interaction, cl, xp, vp, r_cut)
            @test length(cl.next_index) == length(xp)
            empty!(cl)
            @test isempty(cl.next_index)
            @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
            @test run_naive.interaction ≈ run_cl.interaction rtol=1e-14  # basically the same result
        end

        @testset "Using foreach_pair" begin
            run_cl = compute_interaction_foreach_pair(f_interaction, cl, xp, vp, r_cut)
            @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
            @test run_naive.interaction ≈ run_cl.interaction rtol=1e-13  # basically the same result
        end

        @testset "Using foreach_pair (batched)" begin
            run_cl = compute_interaction_foreach_pair(f_interaction, cl, xp, vp, r_cut; batch_size = Val(4))
            @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
            @test run_naive.interaction ≈ run_cl.interaction rtol=1e-13  # basically the same result
        end

        @testset "Using foreach_pair (folded)" begin
            # This assumes points are in [0, L]
            run_cl = compute_interaction_foreach_pair(f_interaction, cl, xp, vp, r_cut; folded = Val(true))
            @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
            @test run_naive.interaction ≈ run_cl.interaction rtol=1e-13  # basically the same result
        end

        @testset "Using foreach_pair (PseudoGPU)" begin
            backend = PseudoGPU()
            cl_gpu = construct_cell_list(r_cut, Ls; backend, nsubdiv = Val(nsubdiv))
            xp_gpu = adapt(backend, xp)
            vp_gpu = adapt(backend, vp)
            run_cl = compute_interaction_foreach_pair(f_interaction, cl_gpu, xp_gpu, vp_gpu, r_cut)
            @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
            @test run_naive.interaction ≈ run_cl.interaction rtol=1e-13  # basically the same result
        end

        @testset "Using foreach_source" begin
            run_cl = compute_interaction_foreach_source(f_interaction, cl, xp, vp, r_cut)
            @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
            @test run_naive.interaction ≈ run_cl.interaction rtol=1e-13  # basically the same result
        end
    end

    # Benchmarks
    # cl = construct_cell_list(r_cut, Ls; nsubdiv = Val(2))
    # print("nearby_elements (serial):\t")
    # @btime compute_interaction_nearby_elements($f_interaction, $cl, $xp, $vp, $r_cut)  # 2.7 ms (serial)
    # print("foreach_pair:\t")
    # @btime compute_interaction_foreach_pair($f_interaction, $cl, $xp, $vp, $r_cut)  # 450μs (4 threads)
    # print("foreach_source:\t")
    # @btime compute_interaction_foreach_source($f_interaction, $cl, $xp, $vp, $r_cut)  # 450μs (4 threads)

    nothing
end

@testset "PeriodicCellList" begin
    test_cell_lists()
end
