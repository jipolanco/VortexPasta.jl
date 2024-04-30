using LinearAlgebra: ×, norm
using SpecialFunctions: erfc, erf
using ..Filaments: deperiodise_separation, Segment, segments

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

ewald_screening_function(::Velocity, ::ShortRange, αr::Real) = erfc(αr) + 2αr / sqrt(π) * exp(-αr^2)
ewald_screening_function(::Velocity, ::ShortRange,   ::Zero) = 1

ewald_screening_function(::Velocity, ::LongRange, αr::Real) = erf(αr) - 2αr / sqrt(π) * exp(-αr^2)
ewald_screening_function(::Velocity, ::LongRange,   ::Zero) = Zero()

ewald_screening_function(::Streamfunction, ::ShortRange, αr::Real) = erfc(αr)
ewald_screening_function(::Streamfunction, ::ShortRange,   ::Zero) = 1

ewald_screening_function(::Streamfunction, ::LongRange, αr::Real) = erf(αr)
ewald_screening_function(::Streamfunction, ::LongRange,   ::Zero) = Zero()

# This simply corresponds to the unmodified BS integrals.
ewald_screening_function(::Velocity, ::FullIntegrand, αr) = 1
ewald_screening_function(::Streamfunction, ::FullIntegrand, αr) = 1

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
    @inbounds for i ∈ eachindex(Xs, vs)
        x⃗ = Xs[i]
        sa = Segment(f, ifelse(i == firstindex(segs), lastindex(segs), i - 1))  # segment i - 1 (with periodic wrapping)
        sb = Segment(f, i)  # segment i
        u⃗a = integrate_biot_savart(quantity, LongRange(), sa, x⃗, params; Lhs, rcut² = nothing)
        u⃗b = integrate_biot_savart(quantity, LongRange(), sb, x⃗, params; Lhs, rcut² = nothing)
        vs[i] = vs[i] - prefactor * (u⃗a + u⃗b)
    end
    vs
end

function remove_long_range_self_interaction!(
        vs::AbstractVector{<:VectorOfVec},
        fs::VectorOfFilaments,
        args...,
    )
    Threads.@threads :static for i ∈ eachindex(vs, fs)
        @inbounds v, f = vs[i], fs[i]
        remove_long_range_self_interaction!(v, f, args...)
    end
    vs
end

function add_short_range_fields!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::ShortRangeCache,
        fs::VectorOfFilaments;
        kws...,
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    # We parallelise at the level of the list of filaments.
    # This is good when there are many filaments.
    # Moreover, the :static scheduling option makes sense when roughly all the filaments
    # have the same number of points. If that's not the case, it may be worth it to either
    # (1) reorder filaments such that every "chunk" has oroughly the same number of points,
    # or (2) use :dynamic scheduling.
    Threads.@threads :static for i ∈ eachindex(fs)
        fields_i = map(us -> us[i], fields)  # velocity/streamfunction of i-th filament
        add_short_range_fields!(fields_i, cache, fs[i]; kws...)
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
    (; quad, lia_segment_fraction,) = params
    (; Γ, a, Δ, Ls, quad_near_singularity,) = params.common
    prefactor = Γ / (4π)
    Lhs = map(L -> L / 2, Ls)

    Xs = nodes(f)
    for us ∈ values(fields)
        eachindex(us) == eachindex(Xs) ||
            throw(DimensionMismatch("vector has wrong length"))
    end

    segs = segments(f)
    nonlia_lims = nonlia_integration_limits(lia_segment_fraction)
    ps = _fields_to_pairs(fields)

    for i ∈ eachindex(Xs)
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

        # Then include the short-range (but non-local) effect of all nearby charges
        # (which are located on the quadrature points of all nearby segments,
        # excluding the local segments `sa` and `sb`).
        for charge ∈ nearby_charges(cache, x⃗)
            s⃗, q⃗, seg = charge
            is_local_segment = seg === sa || seg === sb
            qs⃗′ = real(q⃗)  # just in case data is complex (case of NaiveShortRangeBackend)
            vecs_i = map(vecs_i) do (quantity, u⃗)
                @inline
                if !is_local_segment
                    u⃗ = u⃗ + biot_savart_contribution(
                        quantity, ShortRange(), params.common, x⃗, s⃗, qs⃗′;
                        Lhs, rcut² = params.rcut_sq,
                    )
                end
                quantity => u⃗
            end
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

include("lia.jl")  # defines local_self_induced_velocity (computation of LIA term)
include("naive.jl")
include("cell_lists.jl")
