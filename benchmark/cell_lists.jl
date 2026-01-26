module CellListsBenchmarks

ENV["POCL_AFFINITY"] = 1  # not sure if this is useful
ENV["POCL_CPU_MAX_CU_COUNT"] = Threads.nthreads()  # limit number of CPU threads used by OpenCLBackend (based on POCL)
ENV["POCL_WORK_GROUP_METHOD"] = "cbs"  # might help avoid crashes (https://github.com/pocl/pocl/issues/1971#issuecomment-3062532073)

using VortexPasta.CellLists
using StaticArrays: SVector
using StableRNGs: StableRNG
using KernelAbstractions: KernelAbstractions as KA, CPU
using SIMD: SIMD
using BenchmarkTools
using Adapt: adapt
using OpenCL, pocl_jll

using ThreadPinning
pinthreads(:cores)
threadinfo()

# Print OpenCL information
OpenCL.versioninfo()
@show cl.platform()
@show cl.device()

@inline function deperiodise_separation_folded(r⃗::Vec, Ls, Ls_half) where {Vec}
    map(deperiodise_separation_folded, r⃗, Vec(Ls), Vec(Ls_half))::typeof(r⃗)
end

@inline function deperiodise_separation_folded(r::Real, L::Real, Lh::Real)
    # @assert -L < r < L  # this is true if both points x and y are in [0, L) (r = x - y)
    r = ifelse(r ≥ +Lh, r - L, r)
    r = ifelse(r < -Lh, r + L, r)
    # @assert abs(r) ≤ Lhalf
    r
end

@inline function deperiodise_separation_folded(r::SIMD.Vec, L::Real, Lh::Real)
    # @assert -L < r < L  # this is true if both points x and y are in [0, L) (r = x - y)
    r = SIMD.vifelse(r ≥ +Lh, r - L, r)
    r = SIMD.vifelse(r < -Lh, r + L, r)
    # @assert abs(r) ≤ Lhalf
    r
end

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
            # Assume both points are already in the main periodic cell.
            r⃗ = deperiodise_separation_folded(x⃗ - y⃗, Ls, Ls_half)
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
        CellLists.foreach_source(cl, x⃗; folded = Val(true)) do j
            @inbounds y⃗ = xp[j]
            # Assume both points are already in the main periodic cell.
            r⃗ = deperiodise_separation_folded(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            @inbounds if r² <= r²_cut
                wp[i] += f(vp[i], vp[j], r⃗, r²)
            end
        end
    end
    wp
end

function interactions_foreach_pair!(f::F, cl::PeriodicCellList, wp, xp, vp, r_cut; simd = false, sort_points = true) where {F}
    (; Ls,) = cl
    Ls_half = Ls ./ 2
    r²_cut = r_cut^2
    fill!(wp, 0)
    if simd
        CellLists.foreach_pair(cl, xp; folded = Val(true), sort_points, batch_size = Val(4)) do x⃗, i, js, m
            # @inbounds x⃗ = xp[i]
            W = length(js)
            j = SIMD.Vec(js)
            N = length(x⃗)  # number of dimensions (= 3)
            # xp_tup = StructArrays.components(xp)::NTuple{N}
            ys = ntuple(Val(N)) do d
                SIMD.Vec(map(l -> @inbounds(xp[l][d]), js))
            end
            # ys = map(x -> @inbounds(x[j]), xp_tup)::NTuple{N, SIMD.Vec}
            rs = ntuple(Val(N)) do d
                @inbounds deperiodise_separation_folded(x⃗[d] - ys[d], Ls[d], Ls_half[d])::SIMD.Vec
            end
            r² = sum(abs2, rs)::SIMD.Vec
            mask = SIMD.Vec(ntuple(identity, Val(W))) ≤ m
            mask = mask & (r² ≤ r²_cut)
            if any(mask)
                vi = @inbounds vp[i]
                vj = @inbounds vp[j]
                ws = f(vi, vj, rs, r²)
                @inbounds wp[i] += sum(SIMD.vifelse(mask, ws, zero(ws)))
            end
        end
    else
        CellLists.foreach_pair(cl, xp; folded = Val(true), sort_points) do x⃗, i, j
            @inbounds y⃗ = xp[j]
            # Assume both points are already in the main periodic cell.
            r⃗ = deperiodise_separation_folded(x⃗ - y⃗, Ls, Ls_half)
            r² = sum(abs2, r⃗)
            @inbounds if r² <= r²_cut
                wp[i] += f(vp[i], vp[j], r⃗, r²)
            end
        end
    end
    KA.synchronize(cl.backend)
    wp
end

function benchmark_set_elements!(cl, xp)
    CellLists.set_elements!(cl, xp; folded = Val(true))
    KA.synchronize(cl.backend)
    nothing
end

function main()
    suite = BenchmarkGroup()

    T = Float64
    Ls = (2π, 3π, 4π)
    r_cut = 0.1
    rs_cut = (r_cut, r_cut, r_cut)

    Np = 400_000
    rng = StableRNG(42)

    # Random locations and values within the domain (already folded into [0, L]^3)
    xp_cpu = rand(rng, SVector{3, T}, Np)
    for i in eachindex(xp_cpu)
        xp_cpu[i] = xp_cpu[i] .* Ls
    end
    vp_cpu = randn(rng, T, Np)
    wp_cpu = similar(vp_cpu)  # this is the "influence" on a point of all neighbouring points

    # Arbitrary interaction function
    function f_interaction(u::T, v, r⃗, r²) where {T}
        r_c = T(2π / 100)  # to avoid singularity
        r = sqrt(r²)
        u * v / (r + r_c) * exp(-r² / r_c^2)
    end

    backends = [
        "CPU" => CPU(),
        "OpenCLBackend" => OpenCLBackend(),
    ]

    for (backend_name, backend) in backends, nsubdiv in 1:2
        cl = PeriodicCellList(backend, rs_cut, Ls, Val(nsubdiv))
        xp = adapt(backend, xp_cpu)
        vp = adapt(backend, vp_cpu)
        wp = adapt(backend, wp_cpu)
        CellLists.set_elements!(cl, xp; folded = Val(true))
        if backend isa CPU
            wp_a = copy(interactions_foreach_pair!(f_interaction, cl, wp, xp, vp, r_cut; sort_points = true))
            wp_b = copy(interactions_foreach_pair!(f_interaction, cl, wp, xp, vp, r_cut; simd = true))
            wp_c = copy(interactions_foreach_pair!(f_interaction, cl, wp, xp, vp, r_cut; sort_points = false))
            @assert wp_a ≈ wp_b  # verification
            @assert wp_a ≈ wp_c  # verification
        end
        sub = suite[backend_name]["nsubdiv = $nsubdiv"]
        sub["set_elements!"] = @benchmarkable benchmark_set_elements!($cl, $xp)
        if backend isa CPU
            sub["iterator_interface"] = @benchmarkable interactions_iterator_interface!($f_interaction, $cl, $wp, $xp, $vp, $r_cut)
            sub["foreach_source"] = @benchmarkable interactions_foreach_source!($f_interaction, $cl, $wp, $xp, $vp, $r_cut)
        end
        sub["foreach_pair (unsorted)"] = @benchmarkable interactions_foreach_pair!($f_interaction, $cl, $wp, $xp, $vp, $r_cut; sort_points = false)
        sub["foreach_pair (sorted)"] = @benchmarkable interactions_foreach_pair!($f_interaction, $cl, $wp, $xp, $vp, $r_cut; sort_points = true)
        if backend isa CPU
            sub["foreach_pair (SIMD/sorted)"] = @benchmarkable interactions_foreach_pair!($f_interaction, $cl, $wp, $xp, $vp, $r_cut; simd = true, sort_points = true)
            sub["foreach_pair (SIMD/unsorted)"] = @benchmarkable interactions_foreach_pair!($f_interaction, $cl, $wp, $xp, $vp, $r_cut; simd = true, sort_points = false)
        end
    end

    suite
end

end
