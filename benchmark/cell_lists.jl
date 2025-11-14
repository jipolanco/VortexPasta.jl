module CellListsBenchmarks

using VortexPasta.CellLists
using StaticArrays: SVector
using StableRNGs: StableRNG
using KernelAbstractions: CPU
using BenchmarkTools

@inline function deperiodise_separation(r⃗::Vec, Ls, Ls_half) where {Vec}
    map(deperiodise_separation, r⃗, Vec(Ls), Vec(Ls_half))::typeof(r⃗)
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

function interactions_iterator_interface!(f::F, cl::PeriodicCellList, wp, xp, vp, r_cut) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    fill!(wp, 0)
    @inbounds Threads.@threads for i in eachindex(xp, vp)
        x⃗ = xp[i]
        for j in CellLists.nearby_elements(cl, x⃗)
            @inbounds y⃗ = xp[j]
            r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            @inbounds if r² <= r²_cut
                wp[i] += f(vp[i], vp[j], r⃗, r²)
            end
        end
    end
    wp
end

function interactions_foreach_source!(f::F, cl::PeriodicCellList, wp, xp, vp, r_cut) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    fill!(wp, 0)
    @inbounds Threads.@threads for i in eachindex(xp, vp)
        x⃗ = xp[i]
        CellLists.foreach_source(cl, x⃗) do j
            @inbounds y⃗ = xp[j]
            r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            @inbounds if r² <= r²_cut
                wp[i] += f(vp[i], vp[j], r⃗, r²)
            end
        end
    end
    wp
end

function interactions_foreach_pair!(f::F, cl::PeriodicCellList, wp, xp, vp, r_cut) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    fill!(wp, 0)
    CellLists.foreach_pair(cl, xp) do x⃗, i, j
        @inbounds y⃗ = xp[j]
        r⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Ls_half)
        r² = sum(abs2, r⃗)
        @inbounds if r² <= r²_cut
            wp[i] += f(vp[i], vp[j], r⃗, r²)
        end
    end
    wp
end

function main()
    suite = BenchmarkGroup()

    T = Float64
    Ls = (2π, 3π, 4π)
    r_cut = 0.1
    rs_cut = (r_cut, r_cut, r_cut)

    Np = 400_000
    rng = StableRNG(42)

    # Random locations and values within the domain
    xp = rand(rng, SVector{3, T}, Np)
    for i in eachindex(xp)
        xp[i] = xp[i] .* Ls
    end
    vp = randn(rng, T, Np)
    wp = similar(vp)  # this is the "influence" on a point of all neighbouring points

    # Arbitrary interaction function
    function f_interaction(u, v, r⃗, r²)
        r_c = oftype(r², 2π / 100)  # to avoid singularity
        r = sqrt(r²)
        u * v / (r + r_c) * exp(-r² / r_c^2)
    end

    for nsubdiv in 1:2
        cl = PeriodicCellList(CPU(), rs_cut, Ls, Val(nsubdiv))
        CellLists.set_elements!(cl, xp)
        sub = suite["nsubdiv = $nsubdiv"]
        sub["set_elements!"] = @benchmarkable CellLists.set_elements!($cl, $xp)
        sub["iterator_interface"] = @benchmarkable interactions_iterator_interface!($f_interaction, $cl, $wp, $xp, $vp, $r_cut)
        sub["foreach_source"] = @benchmarkable interactions_foreach_source!($f_interaction, $cl, $wp, $xp, $vp, $r_cut)
        sub["foreach_pair"] = @benchmarkable interactions_foreach_pair!($f_interaction, $cl, $wp, $xp, $vp, $r_cut)
    end


    suite
end

end
