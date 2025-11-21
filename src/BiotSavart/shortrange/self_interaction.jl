# TODO: these integrals could be directly computed from PointData (e.g. on GPU), as we don't
# need curvatures or anything complicated.

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

## ========================================================================================== ##
## Computation of local integral when lia_segment_fraction < 1

function add_local_integrals!(
        vs::AbstractVector{<:VectorOfVec}, fs::VectorOfFilaments, quantity::OutputField, params::ParamsCommon;
        lia_segment_fraction
    )
    lia_segment_fraction === nothing && return vs
    chunks = FilamentChunkIterator(fs)
    @sync for chunk in chunks
        isempty(chunk) && continue  # don't spawn a task if it will do no work
        Threads.@spawn for i in chunk
            @inbounds v, f = vs[i], fs[i]
            add_local_integrals!(v, f, quantity, params; lia_segment_fraction)
        end
    end
    vs
end

function add_local_integrals!(
        vs::VectorOfVec, f::ClosedFilament, quantity::OutputField, params::ParamsCommon;
        lia_segment_fraction
    )
    lia_segment_fraction === nothing && return vs
    (; quad_near_singularity) = params
    lims = nonlia_integration_limits(lia_segment_fraction)
    Xs = nodes(f)
    segs = segments(f)
    Lhs = map(L -> L / 2, params.Ls)
    prefactor = params.Γ / (4π)
    @inbounds for i in eachindex(Xs, vs)
        x⃗ = Xs[i]
        sa = Segment(f, ifelse(i == firstindex(segs), lastindex(segs), i - 1))  # segment i - 1 (with periodic wrapping)
        sb = Segment(f, i)  # segment i
        u⃗a = integrate_biot_savart(quantity, FullIntegrand(), sa, x⃗, params; Lhs, limits = lims[1], quad = quad_near_singularity)
        u⃗b = integrate_biot_savart(quantity, FullIntegrand(), sb, x⃗, params; Lhs, limits = lims[2], quad = quad_near_singularity)
        vs[i] = vs[i] + prefactor * (u⃗a + u⃗b)
    end
    vs
end

function add_local_integrals!(
        fields::NamedTuple, params::ParamsBiotSavart, fs::VectorOfFilaments,
    )
    (; lia_segment_fraction,) = params.shortrange
    lia_segment_fraction === nothing && return fields  # nothing to do
    ps = _fields_to_pairs(fields)
    foreach(ps) do (quantity, vs)
        add_local_integrals!(vs, fs, quantity, params.common; lia_segment_fraction)
    end
    fields
end
