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
process_point_charges!(::ShortRangeCache, ::PointData) = nothing  # can be overriden by the backend

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
        params::ParamsShortRange,
        x⃗::Vec3,                    # point where quantity is computed
        s⃗::Vec3,                    # curve location s⃗(t)
        qs⃗′::Vec3{<:Real};          # curve derivative s⃗′(t) (optionally multiplied by a quadrature weight)
        Lhs = map(L -> L / 2, params.common.Ls),  # this allows to precompute Ls / 2
    )
    (; common, rcut_sq,) = params
    (; Ls, α,) = common
    r⃗ = deperiodise_separation(x⃗ - s⃗, Ls, Lhs)
    r² = sum(abs2, r⃗)
    if component === ShortRange() && r² > rcut_sq
        # The short-range backend may allow some pair interactions beyond the cutoff.
        return zero(r⃗)
    end
    @assert typeof(r⃗) === typeof(qs⃗′)  # should be true for type stability reasons...
    modified_bs_integrand(quantity, component, qs⃗′, r⃗, r², α)
end

# Compute Biot-Savart integral over a single filament segment.
# Note: this doesn't include the prefactor Γ/4π.
function integrate_biot_savart(
        quantity::OutputField,
        component::EwaldComponent,
        seg::Segment,
        x⃗::Vec3,
        params::ParamsShortRange;
        Lhs = map(L -> L / 2, params.common.Ls),  # this allows to precompute Ls / 2
        limits = nothing,
    )
    integrate(seg, params.quad; limits) do seg, ζ
        s⃗ = seg(ζ)
        s⃗′ = seg(ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
        biot_savart_contribution(quantity, component, params, x⃗, s⃗, s⃗′; Lhs)
    end
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
        f::AbstractFilament;
        LIA::Val{_LIA} = Val(true),  # can be used to disable LIA
    ) where {Names, N, V <: VectorOfVec, _LIA}
    ps = _fields_to_pairs(fields)

    (; params,) = cache
    (; quad, regularise_binormal, lia_segment_fraction,) = params
    (; Γ, a, Δ, Ls,) = params.common
    prefactor = Γ / (4π)
    Lhs = map(L -> L / 2, Ls)

    Xs = nodes(f)
    for (_, us) ∈ ps
        eachindex(us) == eachindex(Xs) || throw(DimensionMismatch(
            "vector has wrong length"
        ))
    end

    segs = segments(f)
    nonlia_lims = nonlia_integration_limits(lia_segment_fraction)

    part = indices_per_thread(eachindex(Xs))

    @sync for inds ∈ part
        Threads.@spawn for i ∈ inds
            x⃗ = Xs[i]

            # Start with the "singular" region (i.e. the segments which include x⃗: `sa` and `sb`).
            # We first subtract the effect of the long-range estimation, and then replace it with the LIA term.
            # Removing the long-range integral is needed to have a total velocity which does not
            # depend on the (unphysical) Ewald parameter α.
            # Note that the integral with the long-range kernel is *not* singular (since it's a
            # smoothing kernel), so there's no problem with evaluating this integral.
            # Note: we use a prefactor of 1 (instead of Γ/4π), since we intend to add the prefactor later.
            # We also subtract short-range contributions, since we compute them further below
            # and we're not supposed to include them in the computations.
            # So, in the end, we use FullIntegrand to add both contributions.
            sa = Segment(f, ifelse(i == firstindex(segs), lastindex(segs), i - 1))  # segment i - 1 (with periodic wrapping)
            sb = Segment(f, i)  # segment i

            vecs_i = map(ps) do (quantity, _)
                u⃗ = -(
                    + integrate_biot_savart(quantity, FullIntegrand(), sa, x⃗, params; Lhs)
                    + integrate_biot_savart(quantity, FullIntegrand(), sb, x⃗, params; Lhs)
                )
                if lia_segment_fraction !== nothing
                    # In this case we need to include the integral over a fraction of the local segments.
                    u⃗ = u⃗ + (
                        + integrate_biot_savart(quantity, FullIntegrand(), sa, x⃗, params; Lhs, limits = nonlia_lims[1])
                        + integrate_biot_savart(quantity, FullIntegrand(), sb, x⃗, params; Lhs, limits = nonlia_lims[2])
                    )
                end
                if _LIA
                    u⃗ = u⃗ + local_self_induced(
                        quantity, f, i, one(prefactor);
                        a, Δ, quad, regularise_binormal,
                        segment_fraction = lia_segment_fraction,
                    )
                end
                quantity => u⃗
            end

            # Then integrate short-range effect of nearby charges.
            # Note that this includes the contributions of the two segments within the singular
            # region which we're not supposed to include (we subtracted the same contributions
            # above, so it cancels out).
            for charge ∈ nearby_charges(cache, x⃗)
                s⃗, q⃗ = charge
                qs⃗′ = real(q⃗)  # just in case data is complex (case of NaiveShortRangeBackend)
                vecs_i = map(vecs_i) do (quantity, u⃗)
                    u⃗ = u⃗ + biot_savart_contribution(quantity, ShortRange(), params, x⃗, s⃗, qs⃗′; Lhs)
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
    end

    fields
end

include("lia.jl")  # defines local_self_induced_velocity (computation of LIA term)
include("naive.jl")
include("cell_lists.jl")
