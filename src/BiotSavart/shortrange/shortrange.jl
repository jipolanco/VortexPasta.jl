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

abstract type EwaldComponent end
struct ShortRange <: EwaldComponent end
struct LongRange <: EwaldComponent end
struct FullIntegrand <: EwaldComponent end  # ShortRange + LongRange

erfc(x::SIMD.Vec) = one(x) - verf(x)
erfc(x::AbstractFloat) = SpecialFunctions.erfc(x)

erf(x::AbstractFloat) = SpecialFunctions.erf(x)
erf(::Zero) = Zero()
erfc(::Zero) = One()

@inline two_over_sqrt_pi(::SIMD.Vec{W,T}) where {W,T} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::Zero) = Zero()  # we don't really care about this value; it gets multiplied by Zero() anyway

ewald_screening_function(::Velocity,       ::ShortRange, αr::Real) = erfc(αr) + two_over_sqrt_pi(αr) * αr * exp(-αr^2)
ewald_screening_function(::Streamfunction, ::ShortRange, αr::Real) = erfc(αr)

ewald_screening_function(::Velocity,       ::LongRange, αr::Real) = erf(αr) - two_over_sqrt_pi(αr) * αr * exp(-αr^2)
ewald_screening_function(::Streamfunction, ::LongRange, αr::Real) = erf(αr)

# This simply corresponds to the unmodified BS integrals.
ewald_screening_function(::Velocity,       ::FullIntegrand, αr) = true
ewald_screening_function(::Streamfunction, ::FullIntegrand, αr) = true

biot_savart_integrand(::Velocity, s⃗′, r⃗, r) = (s⃗′ × r⃗) / r^3
biot_savart_integrand(::Streamfunction, s⃗′, r⃗, r) = s⃗′ / r

# Evaluation of long-range integrands at r = 0.
# These are non-singular and can be obtained from Taylor expansions close to r = 0.
# Currently we don't use them but we keep them here just in case.
long_range_bs_integrand_at_zero(::Velocity, α, s⃗′) = zero(s⃗′)
long_range_bs_integrand_at_zero(::Streamfunction, α, s⃗′) = 2 * α / sqrt(π) * s⃗′

@inline function modified_bs_integrand(quantity, component, s⃗′, r⃗, r², α)
    # This is generally not reached as we don't evaluate exactly on nodes.
    # But it could be useful in the future.
    if component === LongRange() && iszero(r²)
        # In this case the result of `biot_savart_integrand` is `NaN`, but the actual
        # integrand when multiplied by the splitting function is well defined.
        return oftype(s⃗′, long_range_bs_integrand_at_zero(quantity, α, s⃗′))
    end
    r = sqrt(r²)
    g = ewald_screening_function(quantity, component, α * r)
    w⃗ = biot_savart_integrand(quantity, s⃗′, r⃗, r)
    oftype(s⃗′, g * w⃗)  # we assume w⃗ is a vector...
end

# Compute the contribution of a single filament point to the Biot-Savart integral on a given
# point x⃗. Note: this doesn't include the prefactor Γ/4π.
@inline function biot_savart_contribution(
        quantity::OutputField,      # Velocity() or Streamfunction()
        component::EwaldComponent,  # ShortRange() or LongRange()
        params::ParamsCommon,
        x⃗::V,    # point where quantity is computed
        s⃗::V,    # curve location s⃗(t)
        qs⃗′::V;  # curve derivative s⃗′(t) (optionally multiplied by a quadrature weight)
        Lhs = map(L -> L / 2, params.Ls),  # this allows to precompute Ls / 2
        rcut² = nothing,
    ) where {V <: Vec3{<:Real}}
    # Cut-off distance should be given if and only if we're computing the short-range
    # component.
    @assert (rcut² !== nothing) === (component === ShortRange())
    (; Ls, α,) = params
    r⃗ = deperiodise_separation(x⃗ - s⃗, Ls, Lhs)
    r² = sum(abs2, r⃗)
    if rcut² !== nothing && r² > rcut²
        # The short-range backend may allow some pair interactions beyond the cutoff.
        return zero(r⃗)
    end
    @assert typeof(r⃗) === typeof(qs⃗′)  # should be true for type stability reasons...
    modified_bs_integrand(quantity, component, qs⃗′, r⃗, r², α) :: V
end

# Compute Biot-Savart integral over a single filament segment.
# Note: this doesn't include the prefactor Γ/4π.
function integrate_biot_savart(
        quantity::OutputField,
        component::EwaldComponent,
        seg::Segment,
        x⃗::Vec3,
        params::ParamsCommon;
        quad = params.quad,
        Lhs = map(L -> L / 2, params.Ls),  # this allows to precompute Ls / 2
        rcut² = nothing,
        limits::Union{Nothing, NTuple{2, Real}} = nothing,
    )
    integrate(seg, quad; limits) do seg, ζ
        @inline
        s⃗ = seg(ζ)
        s⃗′ = seg(ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
        biot_savart_contribution(quantity, component, params, x⃗, s⃗, s⃗′; Lhs, rcut²)
    end :: typeof(x⃗)
end

# ==================================================================================================== #

# Try to distribute filaments over different threads so that each thread has approximately
# the same number of filament nodes (discrete points).
# Actually, the number of nodes per filament may be very unequal in practical situations, with
# e.g. a single filament having a lot of points and many other small filaments.
# In this case, we may end up with empty chunks (and thus "inactive" threads), which is ok.
# In fact, since we also use threading across filament nodes (in add_short_range_fields!),
# these threads will also perform work.
struct FilamentChunkIterator{Filaments <: AbstractVector{<:AbstractVector}}
    fs::Filaments
    nchunks::Int
end

FilamentChunkIterator(fs::VectorOfFilaments; nchunks = Threads.nthreads()) =
    FilamentChunkIterator(fs, nchunks)

Base.IteratorSize(::Type{<:FilamentChunkIterator}) = Base.HasLength()
Base.IteratorEltype(::Type{<:FilamentChunkIterator}) = Base.HasEltype()
Base.length(it::FilamentChunkIterator) = it.nchunks
Base.eltype(::FilamentChunkIterator) = typeof(1:10)  # = UnitRange{Int}

function Base.iterate(it::FilamentChunkIterator)
    (; fs,) = it
    Np_total = sum(length, fs)  # total number of filament nodes
    Np_accumulated = zero(Np_total)
    nchunk = 0  # index of current chunk
    j_prev = firstindex(fs) - 1  # last filament of previous chunk
    state = (; nchunk, j_prev, Np_total, Np_accumulated)
    iterate(it, state)
end

function Base.iterate(it::FilamentChunkIterator, state)
    (; fs, nchunks,) = it
    (; nchunk, j_prev, Np_total, Np_accumulated,) = state
    if nchunk == it.nchunks
        return nothing  # we're done iterating
    end
    nchunk += 1
    Np_accumulated_wanted = (nchunk * Np_total) ÷ nchunks  # how many nodes do we want up to this chunk included
    j = j_prev  # where this chunk ends
    i = j + 1   # where this chunk starts
    @inbounds while Np_accumulated < Np_accumulated_wanted && j < lastindex(fs)
        j += 1
        Np_accumulated += length(fs[j])
    end
    j_prev = j
    state_new = (; nchunk, j_prev, Np_total, Np_accumulated)::typeof(state)
    i:j, state_new
end

# ==================================================================================================== #

function add_short_range_fields!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::ShortRangeCache,
        fs::VectorOfFilaments;
        kws...,
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    chunks = FilamentChunkIterator(fs)
    @sync for chunk in chunks
        isempty(chunk) && continue  # don't spawn a task if it will do no work
        Threads.@spawn for i in chunk
            fields_i = map(us -> us[i], fields)  # velocity/streamfunction of i-th filament
            add_short_range_fields!(fields_i, cache, fs[i]; kws...)
        end
    end
    fields
end

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
    (; nodes,) = pointdata
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

@inline function _simd_load_batch(points_t::NTuple{N}, charges_t::NTuple{N}, x⃗, js::NTuple{W}, m::Integer, Ls, Lhs, rcut²) where {N, W}
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
    r⃗s_simd = map(s⃗_vec, Tuple(x⃗), Ls, Lhs) do svec, x, L, Lh
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

function _add_pair_interactions_simd!(
        ::Val{include_local_integration}, outputs, cache
    ) where {include_local_integration}
    (; pointdata, params) = cache
    (; points, charges, node_idx_prev,) = pointdata
    (; quad,) = params
    (; Γ, α, Ls) = params.common
    T = typeof(Γ)
    rcut² = params.rcut_sq
    prefactor = Γ / T(4π)
    Lhs = map(L -> L / 2, Ls)

    W = dynamic(pick_vector_width(T))  # how many simultaneous elements to compute (optimal depends on current CPU)

    points_t = StructArrays.components(points)::NTuple
    charges_t = StructArrays.components(charges)::NTuple

    # We assume both source and destination points have already been folded into the main periodic cell.
    # (This is done in add_point_charges!)
    foreach_pair(cache; batch_size = Val(W), folded = Val(true)) do x⃗, i, js, m
        @inline
        (; mask_simd, r²s_simd, q⃗s_simd, r⃗s_simd) = _simd_load_batch(points_t, charges_t, x⃗, js, m, Ls, Lhs, rcut²)

        mask_final = if include_local_integration
            mask_simd
        else
            # Check whether quadrature point `j` is in one of the two adjacent segments to node `i`.
            Nq = length(quad)
            js_node = (SIMD.Vec(js) + Nq - 1) ÷ Nq  # index of filament nodes right before quadrature points js
            iprev = @inbounds node_idx_prev[i]
            mask1 = mask_simd & (js_node != i)  # exclude quadrature points on the segment to the right of node `i`
            mask2 = mask1 & (js_node != iprev)  # exclude quadrature points on the segment to the left of node `i`
            mask2
        end

        # There's nothing interesting to compute if mask is fully zero.
        if !iszero(SIMD.bitmask(mask_final))
            # The next operations should all take advantage of SIMD.
            rs = sqrt(r²s_simd)
            rs_inv = inv(rs)
            r³s_inv = inv(r²s_simd * rs)
            αr = α * rs
            erfc_αr = erfc(αr)
            αr_sq = αr * αr
            exp_term = two_over_sqrt_pi(αr) * αr * exp(-αr_sq)  # SIMD exp currently doesn't work with CUDA -- `LLVM error: Undefined external symbol "exp"`

            args = (erfc_αr, exp_term, rs_inv, r³s_inv, q⃗s_simd, r⃗s_simd)

            # Note: we can safely sum at index `i` without atomics, since the implementation
            # ensures that only one thread has that index.
            if hasproperty(outputs, :streamfunction)
                δu⃗_simd = short_range_integrand(Streamfunction(), args...)::NTuple{3, SIMD.Vec}
                δu⃗_data = map(δu⃗_simd) do component
                    @inline
                    component = SIMD.vifelse(mask_final, component, zero(component))  # remove contributions of masked elements
                    sum(component)  # reduction operation: sum the W elements
                end
                @inbounds outputs.streamfunction[i] += prefactor * Vec3(δu⃗_data)
            end

            if hasproperty(outputs, :velocity)
                δu⃗_simd = short_range_integrand(Velocity(), args...)::NTuple{3, SIMD.Vec}
                δu⃗_data = map(δu⃗_simd) do component
                    @inline
                    component = SIMD.vifelse(mask_final, component, zero(component))  # remove contributions of masked elements
                    sum(component)  # reduction operation: sum the W elements
                end
                @inbounds outputs.velocity[i] += prefactor * Vec3(δu⃗_data)
            end
        end
    end

    outputs
end

function _add_pair_interactions_nosimd!(
        ::Val{include_local_integration}, outputs, cache
    ) where {include_local_integration}
    (; pointdata, params) = cache
    (; points, charges, node_idx_prev) = pointdata
    (; quad,) = params
    (; Γ, α, Ls) = params.common
    T = typeof(Γ)
    rcut² = params.rcut_sq
    prefactor = Γ / T(4π)
    Lhs = map(L -> L / 2, Ls)

    # We assume both source and destination points have already been folded into the main periodic cell.
    # (This is done in add_point_charges!)
    foreach_pair(cache; folded = Val(true)) do x⃗, i, j
        @inline
        if !include_local_integration
            # Check whether quadrature point `j` is in one of the two adjacent segments to node `i`.
            Nq = length(quad)
            j_node = (j + Nq - 1) ÷ Nq  # index of filament node right before quadrature point j
            if j_node == i
                return  # quadrature point is on the segment to the right of node `i`
            end
            if j_node == @inbounds node_idx_prev[i]
                return  # quadrature point is on the segment to the left of node `i`
            end
        end
        s⃗ = @inbounds points[j]
        qs⃗′ = @inbounds charges[j]
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
            r = sqrt(r²)
            r_inv = 1 / r
            r³_inv = 1 / (r² * r)
            αr = α * r
            erfc_αr = erfc(αr)
            αr_sq = αr * αr
            exp_term = two_over_sqrt_pi(αr) * αr * exp(-αr_sq)
            args = (erfc_αr, exp_term, r_inv, r³_inv, Tuple(qs⃗′), Tuple(r⃗))
            # Note: we can safely sum at index `i` without atomics, since the implementation
            # ensures that only one thread has that index.
            if hasproperty(outputs, :streamfunction)
                @inbounds outputs.streamfunction[i] += prefactor * Vec3(short_range_integrand(Streamfunction(), args...))
            end
            if hasproperty(outputs, :velocity)
                @inbounds outputs.velocity[i] += prefactor * Vec3(short_range_integrand(Velocity(), args...))
            end
        end
    end

    outputs
end

# Note: all input values may be SIMD types (Vec or tuple of Vec).
@inline function short_range_integrand(::Velocity, erfc_αr, exp_term, r_inv::V, r³_inv::V, qs⃗′::NTuple{3, V}, r⃗::NTuple{3, V}) where {V}
    factor = (erfc_αr + exp_term) * r³_inv
    vec = crossprod(qs⃗′, r⃗)
    map(vec) do component
        @inline
        factor * component
    end
end

# Note: V may be a real or a SIMD vector.
@inline function crossprod(u::T, v::T) where {T <: NTuple{3}}
    @inbounds (
        u[2] * v[3] - u[3] * v[2],
        u[3] * v[1] - u[1] * v[3],
        u[1] * v[2] - u[2] * v[1],
    ) :: T
end

@inline function short_range_integrand(::Streamfunction, erfc_αr, exp_term, r_inv, r³_inv, qs⃗′, r⃗)
    factor = erfc_αr * r_inv
    vec = qs⃗′
    map(vec) do component
        @inline
        factor * component
    end
end

include("lia.jl")  # defines local_self_induced_velocity (computation of LIA term)
include("self_interaction.jl")
include("backends/naive.jl")
include("backends/cell_lists.jl")
