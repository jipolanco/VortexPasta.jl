using LinearAlgebra: ×, norm
using StaticArrays: MVector, MMatrix, SVector, SMatrix
using ..Filaments: deperiodise_separation, Segment, segments

# SIMD-specific stuff, in particular for erf implementation
using HostCPUFeatures: pick_vector_width
using Static: dynamic
using SpecialFunctions: SpecialFunctions
using SIMD: SIMD

include("cache_common.jl")

include("SIMDFunctions.jl")
using .SIMDFunctions: verf

"""
    nearby_charges(c::ShortRangeCache, x⃗::Vec3)

Return an iterator over the indices of points that are "close" to the location `x⃗`.
"""
function nearby_charges end

"""
    foreach_charge(f::Function, c::ShortRangeCache, x⃗::Vec3)

Apply function `f` to all charges which are close to point `x⃗`.

The function should be of the form `f(j)` where `j` is the index of a point which is "close"
to `x⃗`.

See also [`CellLists.foreach_source`](@ref), which is used when [`CellListsBackend`](@ref) has been chosen.
"""
function foreach_charge end

"""
    foreach_pair(f::Function, c::ShortRangeCache)

Apply function `f` to all point pairs within the chosen cut-off distance.

The function should be of the form `f(x⃗, i, j)` where `x⃗ = c.pointdata.nodes[i]` is a
destination point, and `(i, j)` are a pair of destination/source indices.

See [`CellLists.foreach_pair`](@ref) for more details (which is called when using the [`CellListsBackend`](@ref)).
"""
function foreach_pair end

function short_range_velocity end

erf(x::SIMD.Vec) = verf(x)
erfc(x::SIMD.Vec) = one(x) - erf(x)

erf(x::AbstractFloat) = SpecialFunctions.erf(x)
erfc(x::AbstractFloat) = SpecialFunctions.erfc(x)

erf(::Zero) = Zero()
erfc(::Zero) = One()

@inline two_over_sqrt_pi(::SIMD.Vec{W, T}) where {W, T} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::Zero) = Zero()  # we don't really care about this value; it gets multiplied by Zero() anyway

# ==================================================================================================== #

# Iterate over a single filament chunk.
# We iterate over filaments fs, possibly starting somewhere in the middle of filament
# fs[begin], and possibly ending somewhere in the middle of filament fs[end].
# Note that fs is a SubArray, i.e. a view of a larger vector of filaments.
struct SingleChunkIterator{Filaments <: AbstractVector{<:AbstractFilament}}
    fs::Filaments         # all filaments (included those not in this chunk)
    inds::UnitRange{Int}  # indices of considered filaments
    a::Int  # index of start node in first filament (inclusive) -> iterate over nodes fs[inds[begin]][a:end]
    b::Int  # index of end node in last filament (inclusive)    -> iterate over nodes fs[inds[end]][begin:b]
end

Base.IteratorSize(::Type{<:SingleChunkIterator}) = Base.HasLength()
Base.length(it::SingleChunkIterator) = length(it.inds)  # = number of filaments included in this chunk
Base.IteratorEltype(::Type{<:FilamentChunkIterator}) = Base.HasEltype()
Base.eltype(::FilamentChunkIterator) = Tuple{Int, UnitRange{Int}, Int}  # = (filament_idx, node_indices, num_nodes_visited)

function Base.iterate(it::SingleChunkIterator)
    (; fs, inds, a, b) = it
    isempty(inds) && return nothing
    # First filament: start from index a
    i = first(inds)
    f = fs[i]
    if length(inds) == 1  # the chunk covers (part of) a single filament only
        node_indices = a:b
    else  # the chunk covers two or more filaments
        node_indices = a:lastindex(f)
    end::UnitRange{Int}
    # Count accumulated number of nodes until the end of the previous filament.
    num_nodes_visited = sum(j -> length(fs[j]), firstindex(fs):(i - 1); init = 0)
    # Include previously visited nodes in this filament (if a > firstindex(f)).
    num_nodes_visited += a - firstindex(f)
    ret = (i, node_indices, num_nodes_visited)
    num_nodes_visited += length(node_indices)  # include nodes visited in this iteration
    state = (firstindex(inds) + 1, num_nodes_visited)  # (index of next filament, number of visited nodes)
    ret, state
end

function Base.iterate(it::SingleChunkIterator, state::Tuple)
    (; fs, inds, b) = it
    j, num_nodes_visited = state
    j == lastindex(inds) + 1 && return nothing  # we're done iterating over filaments
    i = inds[j]
    f = fs[i]
    if j == lastindex(inds)
        node_indices = firstindex(f):b
    else
        node_indices = UnitRange(eachindex(f))
    end::UnitRange{Int}
    ret = (i, node_indices, num_nodes_visited)
    num_nodes_visited += length(node_indices)
    state = (j + 1, num_nodes_visited)
    ret, state
end

# ==================================================================================================== #

# Try to distribute filaments nodes over different threads so that each thread has approximately
# the same number of filament nodes (discrete points). In fact, the number of nodes per
# filament may be very unequal in practical situations, with e.g. a single filament having a
# lot of points and many other small filaments, so this can help with load balancing.
struct FilamentChunkIterator{Filaments <: AbstractVector{<:AbstractVector}}
    fs::Filaments
    nchunks::Int
end

FilamentChunkIterator(fs::VectorOfFilaments; nchunks = Threads.nthreads()) =
    FilamentChunkIterator(fs, nchunks)

Base.IteratorSize(::Type{<:FilamentChunkIterator}) = Base.SizeUnknown()  # the iterator may generate less than `nchunks` elements
Base.IteratorEltype(::Type{<:FilamentChunkIterator}) = Base.HasEltype()
Base.eltype(::FilamentChunkIterator) = typeof((1:10, (2, 4)))  # = Tuple{UnitRange{Int}, Tuple{Int, Int}} = (i:j, (i_node_idx, j_node_idx))

function Base.iterate(it::FilamentChunkIterator)
    (; fs,) = it
    isempty(fs) && return nothing
    Np_total = sum(length, fs)  # total number of filament nodes
    Np_accumulated = zero(Np_total)
    nchunk = 0  # index of current chunk
    Base.require_one_based_indexing(fs)
    Base.require_one_based_indexing(first(fs))
    i_next = firstindex(fs)                  # first filament of next chunk (i)
    i_node_idx_next = firstindex(first(fs))  # first node of filament i_next to be considered (i_node_idx)
    state = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)
    iterate(it, state)
end

function Base.iterate(it::FilamentChunkIterator, state)
    (; fs, nchunks,) = it
    (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated,) = state
    if nchunk == nchunks || i_next == lastindex(fs) + 1
        return nothing  # we're done iterating
    end
    checkbounds(fs[i_next], i_node_idx_next)  # i_node_idx_next is a node index of the next filament
    nchunk += 1
    # Make sure Np_accumulated_wanted is a multiple of a small power of 2. This might help
    # with false sharing issues (not sure).
    p = 16
    Np_accumulated_wanted = ((Np_total * nchunk) ÷ (nchunks * p)) * p
    if nchunk == nchunks
        Np_accumulated_wanted = Np_total  # make sure we iterate over all points by the last iteration
    end
    if Np_accumulated_wanted == Np_accumulated
        # Nothing to do, iterate recursively with new value of nchunk
        state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
        return iterate(it, state_new)
    end
    i = i_next
    i_node_idx = i_node_idx_next
    j = i  # filament where this chunk ends (to be adjusted below)
    j_node_idx = i_node_idx - 1  # node index where this chunk ends (to be adjusted below)
    while Np_accumulated < Np_accumulated_wanted && j ≤ lastindex(fs)
        Np_wanted = Np_accumulated_wanted - Np_accumulated
        # Available nodes in current filament
        Np_available = lastindex(fs[j]) - j_node_idx
        if Np_available > Np_wanted
            # Stop in the middle of the filament
            j_node_idx = j_node_idx + Np_wanted
            Np_accumulated += Np_wanted
            @assert Np_accumulated == Np_accumulated_wanted
            i_next = j  # continue on the same filament
            i_node_idx_next = j_node_idx + 1
            state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
            ret = SingleChunkIterator(fs, i:j, i_node_idx, j_node_idx)
            return ret, state_new
        elseif Np_available == Np_wanted
            # Stop at the end of this filament
            j_node_idx = j_node_idx + Np_wanted
            Np_accumulated += Np_wanted
            @assert Np_accumulated == Np_accumulated_wanted
            i_next = j + 1
            i_node_idx_next = firstindex(fs[j])
            state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
            ret = SingleChunkIterator(fs, i:j, i_node_idx, j_node_idx)
            return ret, state_new
        else
            # Continue iterating over filaments
            Np_accumulated += Np_available
            j += 1  # jump to next filament
            j_node_idx = 0
        end
    end
    @assert Np_accumulated == Np_accumulated_wanted
    @assert j == lastindex(fs)
    j_node_idx = lastindex(fs[j])  # last node of last filament
    i_next = j + 1
    i_node_idx_next = firstindex(fs[j])
    state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
    ret = SingleChunkIterator(fs, i:j, i_node_idx, j_node_idx)
    return ret, state_new
end

# ==================================================================================================== #

"""
    add_pair_interactions!(outputs::NamedTuple, cache::ShortRangeCache)

Compute short-range Biot-Savart interactions between pairs of points.

This function evaluates the influence of nearby quadrature points on filament discretisation
points. Results are _added_ to the original values in `outputs`, which means that one should
set to zero all values _before_ calling this function (unless one wants to keep previous
data, e.g. the velocity associated from other terms not computed here).

The `outputs` argument is expected to have `:velocity` and/or `:streamfunction` fields,
containing linear vectors of the same length as the number of output points nodes
(`cache.pointdata.nodes`).

This function does _not_ compute local interactions, i.e. the term associated to local
curvature effects in the case of velocity.
"""
function add_pair_interactions!(outputs::NamedTuple, cache::ShortRangeCache)
    if cache.params.common.avoid_explicit_erf
        _add_pair_interactions!(Val(true), outputs, cache)
    else
        _add_pair_interactions!(Val(false), outputs, cache)
    end
end

function _add_pair_interactions!(inc::Val{include_local_integration}, outputs, cache) where {include_local_integration}
    (; pointdata, params) = cache
    (; nodes) = pointdata
    (; avoid_explicit_erf) = params.common
    @assert include_local_integration === avoid_explicit_erf  # we include integration over local segment

    @assert 1 ≤ length(outputs) ≤ 2  # velocity and/or streamfunction
    foreach(outputs) do vs
        eachindex(vs) == eachindex(nodes) || throw(ArgumentError("wrong length of output vector"))
    end

    if cache.params.use_simd
        _add_pair_interactions_simd!(inc, outputs, cache)
    else
        _add_pair_interactions_nosimd!(inc, outputs, cache)
    end

    outputs
end

@inline function _simd_deperiodise_separation_folded(r::SIMD.Vec, L::Real, Lh::Real)
    # @assert -L < r < L  # this is true if both points x and y are in [0, L) (r = x - y)
    r = SIMD.vifelse(r ≥ +Lh, r - L, r)
    r = SIMD.vifelse(r < -Lh, r + L, r)
    # @assert abs(r) ≤ Lhalf
    r
end

# Case of non-periodic domains.
@inline _simd_deperiodise_separation_folded(r::SIMD.Vec, L::Infinity, Lh::Infinity) = r

@inline function _simd_load_batch(points_t::NTuple{N}, charges_t::NTuple{N}, x⃗::NTuple{N}, js::NTuple{W}, m::Integer, Ls, Lhs, rcut²) where {N, W}
    # @assert m <= W
    # Note: all indices in `js` are all expected to be valid (even when m < W), so they can
    # be used to index `points` and `charges`.

    Vec = SIMD.Vec
    js_vec = Vec(js)

    # Load W points. Note that, even when m < W, all W indices in `js` are expected to be valid.
    # This is guaranteed by the CellLists.foreach_source implementation in particular.
    s⃗_vec = map(points_t) do xs
        @inline
        @inbounds SIMD.vgather(xs, js_vec)  # load non-contiguous values
    end::NTuple{N, Vec{W}}

    # Check that all points have already been folded in [0, L]
    # (This is assumed further below.)
    # foreach(s⃗_vec, Ls) do xs, L
    #     @assert all(0 ≤ xs) && all(xs < L)
    # end

    # Determine distances between x⃗ and source points s⃗_vec.
    # Note that we want the _minimal_ distance in the periodic lattice.
    r⃗s_simd = map(s⃗_vec, x⃗, Ls, Lhs) do svec, x, L, Lh
        @inline
        # Assuming both source and destination points are in [0, L], we need max two
        # operations to obtain their minimal distance in the periodic lattice.
        local rs = _simd_deperiodise_separation_folded(x - svec, L, Lh)
        # @assert all(-Lh ≤ rs) && all(rs < Lh)
        rs
    end::NTuple{N, Vec{W}}

    q⃗s_simd = map(charges_t) do qs
        @inline
        @inbounds SIMD.vgather(qs, js_vec)
    end::NTuple{N, Vec{W}}

    r²s_simd = sum(abs2, r⃗s_simd)::Vec{W}

    # Mask out elements satisfying at least one of the criteria:
    # - their distance is beyond rcut
    # - their index in 1:W is in (m + 1):W
    one_to_W = Vec(ntuple(identity, Val(W)))
    mask_simd = (r²s_simd ≤ rcut²) & (one_to_W ≤ m)

    (; mask_simd, r²s_simd, q⃗s_simd, r⃗s_simd)
end

# Check if quadrature point j_quad is on the segment to the _right_ of discretisation node `i`.
# Note that j_quad can be a SIMD.Vec.
@inline function _is_on_right_segment(i_node, j_quad, quad::StaticSizeQuadrature)
    Nq = length(quad)
    (Nq * (i_node - 1) < j_quad) & (j_quad ≤ Nq * i_node)
end

function _add_pair_interactions_simd!(
        ::Val{include_local_integration}, outputs::NamedTuple{Names}, cache
    ) where {include_local_integration, Names}
    (; pointdata, params) = cache
    (; points, charges, node_idx_prev,) = pointdata
    (; quad,) = params
    (; Γ, α, Ls) = params.common
    T = typeof(Γ)
    rcut² = params.rcut_sq
    prefactor = Γ / T(4π)
    Lhs = map(L -> L / 2, Ls)
    quantities = NamedTuple{Names}(possible_output_fields())  # e.g. (velocity = Velocity(),)

    W = dynamic(pick_vector_width(T))  # how many simultaneous elements to compute (optimal depends on current CPU)

    points_t = StructArrays.components(points)::NTuple
    charges_t = StructArrays.components(charges)::NTuple

    # We assume both source and destination points have already been folded into the main periodic cell.
    # (This is done in add_point_charges!)
    foreach_pair(cache; batch_size = Val(W), folded = Val(true)) do x⃗, i, js, m
        @inline
        (; mask_simd, r²s_simd, q⃗s_simd, r⃗s_simd) = _simd_load_batch(points_t, charges_t, Tuple(x⃗), js, m, Ls, Lhs, rcut²)

        mask_final = if include_local_integration
            mask_simd
        else
            # Check whether quadrature point `j` is in one of the two adjacent segments to node `i`.
            js_vec = SIMD.Vec(js)
            iprev = @inbounds node_idx_prev[i]
            islocal_r = _is_on_right_segment(i, js_vec, quad)      # quadrature point is on the segment to the right of node `i`
            islocal_l = _is_on_right_segment(iprev, js_vec, quad)  # quadrature point is on the segment to the left of node `i`
            mask_simd & !(islocal_r | islocal_l)
        end

        # There's nothing interesting to compute if mask is fully zero.
        if !iszero(SIMD.bitmask(mask_final))
            # The next operations should all take advantage of SIMD.
            rs = sqrt(r²s_simd)
            rs_inv = inv(rs)
            αr = α * rs
            erfc_αr = erfc(αr)
            exp_term = two_over_sqrt_pi(αr) * αr * exp(-(αr * αr))  # SIMD exp currently doesn't work with CUDA -- `LLVM error: Undefined external symbol "exp"`
            args = (erfc_αr, exp_term, rs_inv, q⃗s_simd, r⃗s_simd)

            foreach(values(outputs), values(quantities)) do vs, quantity
                @inline
                δu⃗_simd = short_range_integrand(quantity, args...)::NTuple{3, SIMD.Vec}
                δu⃗_data = map(δu⃗_simd) do component
                    @inline
                    component = SIMD.vifelse(mask_final, component, zero(component))  # remove contributions of masked elements
                    sum(component)  # reduction operation: sum the W elements
                end
                # Note: we can safely sum at index `i` without atomics, since the implementation
                # ensures that only one thread has that index.
                @inbounds vs[i] += prefactor * Vec3(δu⃗_data)
            end
        end
    end

    outputs
end

function _add_pair_interactions_nosimd!(
        ::Val{include_local_integration}, outputs::NamedTuple{Names}, cache
    ) where {include_local_integration, Names}
    (; pointdata, params) = cache
    (; points, charges, node_idx_prev) = pointdata
    (; quad,) = params
    (; Γ, α, Ls) = params.common
    T = typeof(Γ)
    rcut² = params.rcut_sq
    prefactor = Γ / T(4π)
    Lhs = map(L -> L / 2, Ls)
    quantities = NamedTuple{Names}(possible_output_fields())  # e.g. (velocity = Velocity(),)

    # We assume both source and destination points have already been folded into the main periodic cell.
    # (This is done in add_point_charges!)
    foreach_pair(cache; folded = Val(true)) do x⃗, i, j
        @inline
        if !include_local_integration
            # Check whether quadrature point `j` is in one of the two adjacent segments to node `i`.
            iprev = @inbounds node_idx_prev[i]
            _is_on_right_segment(i, j, quad) && return      # quadrature point is on the segment to the right of node `i`
            _is_on_right_segment(iprev, j, quad) && return  # quadrature point is on the segment to the left of node `i`
        end
        s⃗ = @inbounds Tuple(points[j])
        qs⃗′ = @inbounds Tuple(charges[j])
        N = length(Ls)
        r⃗ = ntuple(Val(N)) do d
            @inline
            local r = @inbounds x⃗[d] - s⃗[d]
            local L, Lh = @inbounds Ls[d], Lhs[d]
            # @assert 0 ≤ x⃗[d] < L
            # @assert 0 ≤ s⃗[d] < L
            # Assuming both source and destination points are in [0, L], we need max two
            # operations to obtain their minimal distance in the periodic lattice.
            r = Filaments.deperiodise_separation_folded(r, L, Lh)
            # @assert -Lh ≤ r < Lh
            r
        end
        r² = sum(abs2, r⃗)
        if r² ≤ rcut²
            assume(r² > 0)  # tell the compiler that we're taking the square root of a positive number
            r = sqrt(r²)
            assume(r > 0)   # tell the compiler that we're not dividing by zero
            r_inv = 1 / r
            αr = α * r
            erfc_αr = erfc(αr)
            exp_term = two_over_sqrt_pi(αr) * αr * exp(-(αr * αr))
            args = (erfc_αr, exp_term, r_inv, qs⃗′, r⃗)
            foreach(values(outputs), values(quantities)) do vs, quantity
                @inline
                # Note: we can safely sum at index `i` without atomics, since the implementation
                # ensures that only one thread has that index.
                @inbounds vs[i] += prefactor * Vec3(short_range_integrand(quantity, args...))
            end
        end
    end

    outputs
end

include("integrands.jl")
include("local_integrals.jl")
include("backends/naive.jl")
include("backends/cell_lists.jl")
