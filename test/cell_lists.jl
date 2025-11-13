using VortexPasta.CellLists
using StaticArrays: SVector
using StableRNGs: StableRNG
using Random
using Test

const Vec3{T} = SVector{3, T}

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

function construct_cell_list(r_cut, Ls; nsubdiv = Val(1))
    rs_cut = map(_ -> r_cut, Ls)
    @inferred PeriodicCellList(rs_cut, Ls, nsubdiv)
end

function compute_interaction_cell_lists(f::F, cl, xp, vp, r_cut) where {F}
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

function test_cell_lists()
    T = Float64
    Ls = (2π, 3π, 4π)
    r_cut = π/4

    Np = 1000
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

    for nsubdiv in (Val(1), Val(2))
        cl = construct_cell_list(r_cut, Ls; nsubdiv)
        @test isempty(cl.next_index)
        run_cl = compute_interaction_cell_lists(f_interaction, cl, xp, vp, r_cut)
        @test length(cl.next_index) == length(xp)
        empty!(cl)
        @test isempty(cl.next_index)
        @test run_naive.n_interactions == run_cl.n_interactions      # same number of considered interactions
        @test run_naive.interaction ≈ run_cl.interaction rtol=1e-15  # basically the same result
    end

    nothing
end

@testset "PeriodicCellList" begin
    test_cell_lists()
end
