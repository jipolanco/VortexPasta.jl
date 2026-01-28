using LinearAlgebra: ×, norm
using StaticArrays: MVector, MMatrix, SVector, SMatrix
using ..Filaments: deperiodise_separation, Segment, segments

# SIMD-specific stuff, in particular for erf implementation
using HostCPUFeatures: pick_vector_width
using Static: dynamic
using SpecialFunctions: SpecialFunctions
using SIMD: SIMD

include("cache_common.jl")

include("SIMDFunctions/SIMDFunctions.jl")
using .SIMDFunctions: SIMDFunctions

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

@inline exp_simd(x::SIMD.Vec) = SIMDFunctions.exp(x)
@inline erf_simd(x::SIMD.Vec) = SIMDFunctions.erf(x)
@inline erfc_simd(x::SIMD.Vec) = one(x) - erf_simd(x)

# Note: even without explicit SIMD, calling SIMD-friendly implementations can enable
# automatic SIMD and thus noticeably improve performance.
@inline exp_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.exp(x)
@inline erf_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.erf(x)
@inline erfc_nosimd(::CPU, x::AbstractFloat) = one(x) - SIMDFunctions.erf(x)

# On GPU we call the functions from Base or SpecialFunctions, since these are usually
# overridden in each GPU implementation (CUDA, ...) and therefore should be fast.
@inline exp_nosimd(::GPU, x::AbstractFloat) = exp(x)
@inline erf_nosimd(::GPU, x::AbstractFloat) = SpecialFunctions.erf(x)
@inline erfc_nosimd(::GPU, x::AbstractFloat) = SpecialFunctions.erfc(x)

@inline exp_nosimd(::KA.Backend, x::Zero) = exp(x)    # = 1 (defined in Constants.jl)
@inline exp_simd(x::Zero) = exp(x)
@inline erf_nosimd(::KA.Backend, ::Zero) = Zero()
@inline erfc_nosimd(::KA.Backend, ::Zero) = One()
@inline erfc_simd(::Zero) = One()

@inline two_over_sqrt_pi(::SIMD.Vec{W, T}) where {W, T} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::Zero) = Zero()  # we don't really care about this value; it gets multiplied by Zero() anyway

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
            erfc_αr = erfc_simd(αr)
            exp_term = two_over_sqrt_pi(αr) * αr * exp_simd(-(αr * αr))
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
    ka_backend = KA.get_backend(cache)
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
            erfc_αr = erfc_nosimd(ka_backend, αr)
            exp_term = two_over_sqrt_pi(αr) * αr * exp_nosimd(ka_backend, -(αr * αr))
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
