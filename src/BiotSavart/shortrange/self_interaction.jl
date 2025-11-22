function remove_self_interaction!(outputs::NamedTuple, cache::ShortRangeCache)
    ka_backend = KA.get_backend(cache)
    if cache.params.common.avoid_explicit_erf
        _remove_self_interaction!(ka_backend, Val(true), outputs, cache)
    else
        _remove_self_interaction!(ka_backend, Val(false), outputs, cache)
    end
    outputs
end

function _remove_self_interaction!(::KA.Backend, avoid_explicit_erf::Val, outputs::NamedTuple{Names}, cache) where {Names}
    (; pointdata, params) = cache
    (; nodes, node_idx_prev, points, charges,) = pointdata
    (; Ls, Γ, α, quad) = params.common

    @assert eachindex(nodes) === eachindex(node_idx_prev)
    foreach(outputs) do vs
        @assert eachindex(vs) === eachindex(nodes)
    end

    T = typeof(Γ)
    prefactor = Γ / T(4π)
    Lhs = map(L -> L / 2, Ls)
    quantities = NamedTuple{Names}(possible_output_fields())  # e.g. (velocity = Velocity(),)

    # This works on CPU (threaded) and GPU.
    AK.foreachindex(nodes; max_tasks = Threads.nthreads(), block_size = 256) do i
        @inline
        x⃗ = @inbounds nodes[i]
        i_prev = @inbounds node_idx_prev[i]
        Nq = length(quad)
        js_right = ntuple(k -> Nq * (i - 1) + k, Val(Nq))      # quadrature points to the right of `i`
        js_left = ntuple(k -> Nq * (i_prev - 1) + k, Val(Nq))  # quadrature points to the left of `i`
        js = (js_left..., js_right...)  # TODO: SIMD over these quadrature points? (optionally?)
        for j in js
            s⃗ = @inbounds Tuple(points[j])
            qs⃗′ = @inbounds Tuple(charges[j])
            N = length(s⃗)
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
            vals = _remove_self_interaction_integral(avoid_explicit_erf, values(quantities), α, r⃗, qs⃗′)::Tuple
            foreach(values(outputs), vals) do vs, v
                @inline
                @inbounds vs[i] -= prefactor * Vec3(v)
            end
        end
    end

    outputs
end

# Compute total BS integral (without Γ/4π prefactor) over local segments.
@inline function _remove_self_interaction_integral(
        avoid_explicit_erf::Val{true}, quantities::Tuple, α, r⃗, qs⃗′,
    )
    r² = sum(abs2, r⃗)
    r = sqrt(r²)
    r_inv = 1 / r
    map(quantities) do quantity
        @inline
        full_integrand(quantity, r_inv, qs⃗′, r⃗)
    end
end

# Compute long-range BS integral (without Γ/4π prefactor) over local segments (needs erf).
@inline function _remove_self_interaction_integral(
        avoid_explicit_erf::Val{false}, quantities::Tuple, α, r⃗, qs⃗′,
    )
    r² = sum(abs2, r⃗)
    r = sqrt(r²)
    r_inv = 1 / r
    αr = α * r
    erf_αr = erf(αr)
    exp_term = two_over_sqrt_pi(αr) * αr * exp(-(αr * αr))
    map(quantities) do quantity
        @inline
        long_range_integrand(quantity, erf_αr, exp_term, r_inv, qs⃗′, r⃗)
    end
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
