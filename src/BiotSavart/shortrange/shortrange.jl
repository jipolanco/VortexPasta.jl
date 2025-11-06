using LinearAlgebra: ×, norm
using StaticArrays: MVector, MMatrix, SVector, SMatrix
using ..Filaments: deperiodise_separation, Segment, segments

# SIMD-specific stuff, in particular for erf implementation
using HostCPUFeatures: pick_vector_width
using Static: dynamic
using SpecialFunctions: SpecialFunctions
using SIMD: SIMD

include("SIMDFunctions.jl")
using .SIMDFunctions: verf

"""
    init_cache_short(
        pc::ParamsCommon, p::ParamsShortRange,
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    ) -> ShortRangeCache

Initialise the cache for the short-range backend defined in `p`.
"""
function init_cache_short end

"""
    process_point_charges!(c::ShortRangeCache, data::PointData)

Process list of point charges.

This is useful for short-range backends like [`CellListsBackend`](@ref), which needs to

Must be called after [`add_point_charges!`](@ref) and before computing any short-range quantities
(using [`add_short_range_fields!`](@ref)).
"""
process_point_charges!(::ShortRangeCache, ::PointData) = nothing  # can be overridden by the backend

"""
    nearby_charges(c::ShortRangeCache, x⃗::Vec3)

Return an iterator over the charges that are "close" to the location `x⃗`.
"""
function nearby_charges end

function short_range_velocity end

abstract type EwaldComponent end
struct ShortRange <: EwaldComponent end
struct LongRange <: EwaldComponent end
struct FullIntegrand <: EwaldComponent end  # ShortRange + LongRange

erfc(x::SIMD.Vec) = one(x) - verf(x)
erfc(x::AbstractFloat) = SpecialFunctions.erfc(x)

erf(x::AbstractFloat) = SpecialFunctions.erf(x)
erf(::Zero) = Zero()
erfc(::Zero) = true  # in the sense of true == 1

two_over_sqrt_pi(::SIMD.Vec{W,T}) where {W,T} = 2 / sqrt(T(π))
two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))
two_over_sqrt_pi(::Zero) = Zero()  # we don't really care about this value; it gets multiplied by Zero() anyway

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
        limits = nothing,
    )
    integrate(seg, quad; limits) do seg, ζ
        @inline
        s⃗ = seg(ζ)
        s⃗′ = seg(ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
        biot_savart_contribution(quantity, component, params, x⃗, s⃗, s⃗′; Lhs, rcut²)
    end :: typeof(x⃗)
end

# This subtracts the local Biot-Savart interaction (the effect of the two adjacent segments
# to each node) which we have included in short and long-range computations, but which
# should be replaced with the local (LIA) term obtained from Taylor expansions to account
# for the integral cut-off near the node.
# This function should be called when the local integrals have _not_ been excluded from
# short-range interactions.
function remove_total_self_interaction!(
        vs::VectorOfVec,
        f::ClosedFilament,
        quantity::OutputField,
        params::ParamsCommon,
    )
    Xs = nodes(f)
    segs = segments(f)
    Lhs = map(L -> L / 2, params.Ls)
    prefactor = params.Γ / (4π)
    # Not sure if we gain much by parallelising here.
    Threads.@threads :dynamic for i in eachindex(Xs, vs)
        @inbounds begin
            x⃗ = Xs[i]
            sa = Segment(f, ifelse(i == firstindex(segs), lastindex(segs), i - 1))  # segment i - 1 (with periodic wrapping)
            sb = Segment(f, i)  # segment i
            u⃗a = integrate_biot_savart(quantity, FullIntegrand(), sa, x⃗, params; Lhs, rcut² = nothing)
            u⃗b = integrate_biot_savart(quantity, FullIntegrand(), sb, x⃗, params; Lhs, rcut² = nothing)
            vs[i] = vs[i] - prefactor * (u⃗a + u⃗b)
        end
    end
    vs
end

# This subtracts the spurious self-interaction term which is implicitly included in
# long-range computations. This corresponds to the integral over the two adjacent segments
# to each node, which should not be included in the total result (since they are replaced by
# the LIA term). This function explicitly computes that integral and subtracts it from the
# result.
#
# This is required because the local term already includes the full (short-range +
# long-range) contribution of these segments. Moreover, long-range
# computations also add the long-range contribution to the
# velocity/streamfunction. Since we don't want to include this contribution
# twice, we subtract it here. Note that without this correction, results will
# depend on the (unphysical) Ewald parameter α. Finally, note that the integral
# with the long-range kernel is *not* singular (since it's a smoothing kernel),
# so there's no problem with evaluating this integral close to x⃗.
function remove_long_range_self_interaction!(
        vs::VectorOfVec,
        f::ClosedFilament,
        quantity::OutputField,
        params::ParamsCommon,
    )
    Xs = nodes(f)
    segs = segments(f)
    Lhs = map(L -> L / 2, params.Ls)
    prefactor = params.Γ / (4π)
    # Note: parallelising here usually leads to worse performance, likely because the work
    # per iteration is small (we might also have false sharing issues due to memory locality).
    @inbounds for i in eachindex(Xs, vs)
        x⃗ = Xs[i]
        sa = Segment(f, ifelse(i == firstindex(segs), lastindex(segs), i - 1))  # segment i - 1 (with periodic wrapping)
        sb = Segment(f, i)  # segment i
        u⃗a = integrate_biot_savart(quantity, LongRange(), sa, x⃗, params; Lhs, rcut² = nothing)
        u⃗b = integrate_biot_savart(quantity, LongRange(), sb, x⃗, params; Lhs, rcut² = nothing)
        vs[i] = vs[i] - prefactor * (u⃗a + u⃗b)
    end
    vs
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

function remove_self_interaction!(
        vs::AbstractVector{<:VectorOfVec},
        fs::VectorOfFilaments, quantity::OutputField, params::ParamsCommon,
    )
    if params.avoid_explicit_erf
        _remove_self_interaction!(Val(true), vs, fs, quantity, params)
    else
        _remove_self_interaction!(Val(false), vs, fs, quantity, params)
    end
end

function _remove_self_interaction!(::Val{total}, vs, fs, quantity, params) where {total}
    chunks = FilamentChunkIterator(fs)
    @sync for chunk in chunks
        isempty(chunk) && continue  # don't spawn a task if it will do no work
        Threads.@spawn for i in chunk
            @inbounds v, f = vs[i], fs[i]
            if total
                # In this case, short-range interactions include the (singular) contribution
                # of the local segments.
                remove_total_self_interaction!(v, f, quantity, params)
            else
                remove_long_range_self_interaction!(v, f, quantity, params)
            end
        end
    end
    vs
end

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
    add_short_range_fields!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::ShortRangeCache,
        f::AbstractFilament;
        LIA = Val(true),
    )

Compute short-range Biot-Savart integrals.

Adds the results onto `fields`.
See [`compute_on_nodes!`](@ref) for more details.

Setting `LIA = Val(false)` allows to disable computation of the localised induction
approximation (LIA) term. In that case, that term should be computed separately using
[`local_self_induced`](@ref).

Before calling this function, one must first call
[`add_point_charges!`](@ref) and then [`process_point_charges!`](@ref).
"""
function add_short_range_fields!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::ShortRangeCache,
        f::ClosedFilament;
        LIA::Val{_LIA} = Val(true),   # can be used to disable LIA
    ) where {Names, N, V <: VectorOfVec, _LIA}
    (; params,) = cache
    if params.common.avoid_explicit_erf
        _add_short_range_fields!(Val(true), fields, cache, f, LIA)
    else
        _add_short_range_fields!(Val(false), fields, cache, f, LIA)
    end
end

function _add_short_range_fields!(::Val{include_local_integration}, fields, cache, f, ::Val{_LIA}) where {include_local_integration, _LIA}
    (; params,) = cache
    (; quad, lia_segment_fraction,) = params
    (; Γ, a, Δ, Ls, quad_near_singularity, avoid_explicit_erf) = params.common
    prefactor = Γ / (4π)
    Lhs = map(L -> L / 2, Ls)
    @assert include_local_integration === avoid_explicit_erf  # we include integration over local segment

    Xs = nodes(f)
    for us ∈ values(fields)
        eachindex(us) == eachindex(Xs) ||
            throw(DimensionMismatch("vector has wrong length"))
    end

    segs = segments(f)
    nonlia_lims = nonlia_integration_limits(lia_segment_fraction)
    ps = _fields_to_pairs(fields)

    # Note: we parallelise here even thought we also parallelised earlier at the level of
    # the filaments. This is ok since we use dynamic threading here (and also at the upper
    # level, since using @spawn also corresponds to dynamic threading)
    Threads.@threads :dynamic for i in eachindex(Xs)
        x⃗ = Xs[i]

        # Determine segments `sa` and `sb` in contact with the singular point x⃗.
        # Then compute the local (LIA) term associated to these segments.
        # If `lia_segment_fraction !== nothing`, we actually compute the local term
        # over a fraction of these segments. Moreover, the rest of the segments
        # (`nonlia_lims`) are taken into account by regular integration using
        # quadratures.
        sa = Segment(f, ifelse(i == firstindex(segs), lastindex(segs), i - 1))  # segment i - 1 (with periodic wrapping)
        sb = Segment(f, i)  # segment i

        vecs_i = map(ps) do (quantity, _)
            @inline
            u⃗ = zero(x⃗)
            if lia_segment_fraction !== nothing
                # In this case we need to include the full BS integral over a fraction of the local segments.
                u⃗ = u⃗ + (
                    integrate_biot_savart(
                        quantity, FullIntegrand(), sa, x⃗, params.common;
                        Lhs, limits = nonlia_lims[1],
                        quad = quad_near_singularity,
                    )
                    +
                    integrate_biot_savart(
                        quantity, FullIntegrand(), sb, x⃗, params.common;
                        Lhs, limits = nonlia_lims[2],
                        quad = quad_near_singularity,
                    )
                )
            end
            if _LIA
                u⃗ = u⃗ + local_self_induced(
                    quantity, f, i, one(prefactor);
                    a, Δ, quad,
                    segment_fraction = lia_segment_fraction,
                )
            end
            quantity => u⃗
        end

        vecs_i = if include_local_integration
            # Then include the short-range effect of all nearby charges, **including the local
            # segments `sa` and `sb`** (which are actually singular, and should be then removed
            # with remove_total_self_interaction).
            add_pair_interactions_shortrange(vecs_i, cache, x⃗, params, nothing, nothing, Lhs)
        else
            # Then include the short-range (but non-local) effect of all nearby charges
            # (which are located on the quadrature points of all nearby segments,
            # excluding the local segments `sa` and `sb`).
            add_pair_interactions_shortrange(vecs_i, cache, x⃗, params, sa, sb, Lhs)
        end

        # Add computed vectors (velocity and/or streamfunction) to corresponding arrays.
        map(ps, vecs_i) do (quantity_p, us), (quantity_i, u⃗)
            @inline
            @assert quantity_p === quantity_i
            us[i] = us[i] + prefactor * u⃗
        end
    end

    fields
end

# Optimised computation of short-range pair interactions.
function add_pair_interactions_shortrange(vecs, cache, x⃗, params::ParamsShortRange, sa, sb, Lhs)
    (; use_simd,) = params
    (; α,) = params.common
    use_simd_actual = use_simd && (α isa AbstractFloat)  # explicit SIMD not used if α == Zero() [non-periodic domains]
    if use_simd_actual
        _add_pair_interactions_shortrange_simd(vecs, cache, x⃗, params, sa, sb, Lhs)
    else
        _add_pair_interactions_shortrange(α, vecs, cache, x⃗, params, sa, sb, Lhs)
    end
end

# Non-periodic case (no Ewald summation, no cutoff distance, no SIMD, naive computation)
function _add_pair_interactions_shortrange(α::Zero, vecs, cache, x⃗, params, sa, sb, Lhs)
    (; Ls,) = params.common
    rcut² = params.rcut_sq
    @assert all(L -> L === Infinity(), Ls)  # infinite non-periodic domain
    @assert all(L -> L === Infinity(), Lhs)
    @assert rcut² === Infinity()
    it = nearby_charges(cache, x⃗)
    for charge ∈ it
        s⃗, q⃗, seg = charge
        if sa !== nothing && sb !== nothing  # if we're excluding local segments from short-range integration
            is_local_segment = seg === sa || seg === sb
            is_local_segment && continue
        end
        qs⃗′ = real(q⃗)  # just in case data is complex
        r⃗ = x⃗ - s⃗
        r² = sum(abs2, r⃗)
        r = sqrt(r²)
        r_inv = 1 / r
        r³_inv = 1 / (r² * r)
        vecs = map(vecs) do (quantity, u⃗)
            δu⃗ = full_integrand(quantity, r_inv, r³_inv, qs⃗′, r⃗)
            u⃗ = u⃗ + δu⃗
            quantity => u⃗
        end
    end
    vecs
end

# Ewald summation without explicit SIMD.
function _add_pair_interactions_shortrange(α::T, vecs, cache, x⃗, params, sa, sb, Lhs) where {T <: AbstractFloat}
    (; Ls,) = params.common
    rcut² = params.rcut_sq
    it = nearby_charges(cache, x⃗)
    for charge ∈ it
        s⃗, q⃗, seg = charge
        is_local_segment = if sa !== nothing && sb !== nothing  # if we're excluding local segments from short-range integration
            seg === sa || seg === sb
        else  # if we want to include local segments in the integration
            false
        end
        qs⃗′ = real(q⃗)  # just in case data is complex
        r⃗ = deperiodise_separation(x⃗ - s⃗, Ls, Lhs)
        r² = sum(abs2, r⃗)
        if is_local_segment || r² > rcut²
            continue
        end
        r = sqrt(r²)
        r_inv = 1 / r
        r³_inv = 1 / (r² * r)
        αr = α * r
        erfc_αr = erfc(αr)
        αr_sq = αr * αr
        exp_term = two_over_sqrt_pi(αr) * αr * exp(-αr_sq)
        vecs = map(vecs) do (quantity, u⃗)
            @inline
            δu⃗ = short_range_integrand(quantity, erfc_αr, exp_term, r_inv, r³_inv, Tuple(qs⃗′), Tuple(r⃗))
            u⃗ = u⃗ + oftype(u⃗, δu⃗)
            quantity => u⃗
        end
    end
    vecs
end

# Non-Ewald (non-periodic) case
full_integrand(::Velocity, r_inv, r³_inv, qs⃗′, r⃗) = r³_inv * (qs⃗′ × r⃗)
full_integrand(::Streamfunction, r_inv, r³_inv, qs⃗′, r⃗) = r_inv * qs⃗′

function mask_to_simd_vec(::Type{SIMD.Vec{W, T}}, mask::Unsigned) where {W, T}
    tup = ntuple(Val(W)) do i
        T((mask >> (i - 1)) & one(mask))
    end
    SIMD.Vec{W, T}(tup)
end

function smatrix_to_simd_vecs(A::SMatrix{W, N, T}) where {W, N, T}
    ntuple(Val(N)) do j
        SIMD.Vec{W, T}(Tuple(A[:, j]))
    end
end

# Periodic case (with Ewald summation)
function _add_pair_interactions_shortrange_simd(vecs, cache, x⃗, params, sa, sb, Lhs)
    (; Ls, α,) = params.common
    @assert !is_open_domain(Ls)
    N = length(x⃗)  # number of dimensions (= 3)
    rcut² = params.rcut_sq
    # with_velocity = any(x -> first(x) === Velocity(), vecs)

    # We use SIMD.jl and MVector/MMatrix to enforce the use of SIMD.
    # This enables important gains in modern CPUs, especially in the computation of erf/erfc.
    T = typeof(α)
    @assert T <: AbstractFloat
    W = dynamic(pick_vector_width(T))  # how many simultaneous elements to compute (optimal depends on current CPU)
    Vec = SIMD.Vec{W, T}

    r²s = MVector{W, T}(undef)

    NW = N * W
    r⃗s = MMatrix{W, N, T, NW}(undef)
    q⃗s = similar(r⃗s)

    it = nearby_charges(cache, x⃗)
    y = iterate(it)

    while y !== nothing
        i = 0
        mask = zero(UInt)  # set to 1 the values that should be computed (0 means ignored)
        @assert 8 * sizeof(mask) ≥ W  # this is basically always the case, but just to be sure...

        # Collect (up to) W charges, to perform operation over W charges at once.
        @inbounds while i < W && y !== nothing
            i += 1
            charge, state = y
            s⃗, q⃗, seg = charge
            is_local_segment = if sa !== nothing && sb !== nothing  # if we're excluding local segments from short-range integration
                seg === sa || seg === sb
            else  # if we want to include local segments in the integration
                false
            end
            r⃗ = deperiodise_separation(x⃗ - s⃗, Ls, Lhs)
            r² = sum(abs2, r⃗)
            for j ∈ eachindex(r⃗, q⃗)
                qj = real(q⃗[j])  # just in case data is complex
                q⃗s[i, j] = qj
                rj = r⃗[j]
                r⃗s[i, j] = rj
            end
            r²s[i] = r²
            # Ignore this element if this is a local segment or if we're beyond the cut-off distance.
            include_element = !is_local_segment && r² ≤ rcut²
            if include_element
                mask = mask | (one(mask) << (i - 1))  # set mask to 1 for this element
            end
            y = iterate(it, state)
        end

        while i < W
            # If we're here, it's because the iterator finished before being able to read W
            # elements. We simply set the missing elements of q⃗s and r⃗s to zero.
            # @assert y === nothing
            i += 1
            @inbounds for j ∈ axes(q⃗s, 2)
                q⃗s[i, j] = 0
                r⃗s[i, j] = 0
            end
        end

        # There's nothing interesting to compute if mask is fully zero.
        iszero(mask) && continue

        # Convert mask to SIMD vector.
        mask_simd = mask_to_simd_vec(SIMD.Vec{W, Bool}, mask)

        # Convert MVector to SIMD vector.
        # The next operations should all take advantage of SIMD.
        r²s_simd = Vec(Tuple(r²s))
        rs = sqrt(r²s_simd)
        rs_inv = inv(rs)
        r³s_inv = inv(r²s_simd * rs)
        αr = α * rs
        erfc_αr = erfc(αr)
        αr_sq = αr * αr
        exp_term = two_over_sqrt_pi(αr) * αr * exp(-αr_sq)

        # Convert SMatrix data onto a tuple of Vec for more explicit SIMD.
        # A single tuple represents a single vector component (e.g. rx, ry and rz).
        # It makes sense to use this to describe vector values (e.g. separation vector,
        # tangent vector, velocity).
        q⃗s_simd = smatrix_to_simd_vecs(SMatrix(q⃗s))
        r⃗s_simd = smatrix_to_simd_vecs(SMatrix(r⃗s))

        vecs = map(vecs) do (quantity, u⃗)
            @inline
            δu⃗_simd = short_range_integrand(quantity, erfc_αr, exp_term, rs_inv, r³s_inv, q⃗s_simd, r⃗s_simd)::NTuple{3, SIMD.Vec}
            δu⃗_data = map(δu⃗_simd) do component
                @inline
                component = SIMD.vifelse(mask_simd, component, zero(component))  # remove contributions of masked elements
                sum(component)  # reduction operation: sum the W elements
            end
            u⃗ = u⃗ + oftype(u⃗, δu⃗_data)
            quantity => u⃗
        end
    end

    vecs
end

# Note: all input values may be SIMD types (Vec or tuple of Vec).
@inline function short_range_integrand(::Velocity, erfc_αr::V, exp_term::V, r_inv::V, r³_inv::V, qs⃗′::NTuple{3, V}, r⃗::NTuple{3, V}) where {V}
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
include("naive.jl")
include("cell_lists.jl")
