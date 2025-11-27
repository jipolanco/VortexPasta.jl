## ========================================================================================== ##
## Remove spurious self-interaction included in other terms

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
        js = (js_left..., js_right...)
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
    assume(r² > 0)  # tell the compiler that we're taking the square root of a positive number
    r = sqrt(r²)
    assume(r > 0)   # tell the compiler that we're not dividing by zero
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
    assume(r² > 0)  # tell the compiler that we're taking the square root of a positive number
    r = sqrt(r²)
    assume(r > 0)  # tell the compiler that we're not dividing by zero
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
## Computation of local term (LIA)

lia_integration_limits(::Nothing) = (nothing, nothing)
lia_integration_limits(γ::Real) = ((one(γ) - γ, one(γ)), (zero(γ), γ))

nonlia_integration_limits(::Nothing) = (nothing, nothing)
nonlia_integration_limits(γ::Real) = ((zero(γ), one(γ) - γ), (γ, one(γ)))

function compute_local_term!(outputs::NamedTuple{Names}, cache::ShortRangeCache) where {Names}
    (; pointdata, params) = cache
    (; derivatives_on_nodes, subsegment_lengths) = pointdata
    (; Γ, a, Δ) = params.common
    @assert eachindex(derivatives_on_nodes[1]) == eachindex(subsegment_lengths[1])

    T = typeof(Γ)
    prefactor = Γ / T(4π)
    quantities = NamedTuple{Names}(possible_output_fields())  # e.g. (velocity = Velocity(),)

    # This works on CPU (threaded) and GPU.
    AK.foreachindex(derivatives_on_nodes[1]; max_tasks = Threads.nthreads(), block_size = 256) do i
        @inline
        # TODO: is it better to precompute sqrt(δ⁻ * δ⁺) on the CPU? (not sure...)
        δ⁻ = @inbounds subsegment_lengths[1][i]
        δ⁺ = @inbounds subsegment_lengths[2][i]
        δ = sqrt(δ⁻ * δ⁺)
        # Derivatives wrt arbitrary parametrisation (=> |s⃗′| ≠ 1)
        s⃗′ = @inbounds derivatives_on_nodes[1][i]  # = Derivative(1)
        s⃗″ = @inbounds derivatives_on_nodes[2][i]  # = Derivative(2)
        s′² = sum(abs2, s⃗′)
        assume(s′² > 0)  # tell the compiler that we're not dividing by zero
        s′_inv = 1 / sqrt(s′²)
        foreach(values(outputs), values(quantities)) do vs, quantity
            @inline
            @inbounds vs[i] = prefactor * _eval_local_term(quantity, δ, s⃗′, s⃗″, s′_inv, a, Δ)::Vec3  # NOTE: we replace old values!
        end
    end

    outputs
end

@inline function _eval_local_term(::Velocity, δ, s⃗′, s⃗″, s′_inv, a, Δ)
    b⃗ = (s⃗′ × s⃗″) * s′_inv^3  # binormal vector (properly normalised) -- see also CurvatureBinormal
    (log(2 * δ / a) - Δ) * b⃗
end

@inline function _eval_local_term(::Streamfunction, δ, s⃗′, s⃗″, s′_inv, a, Δ)
    t̂ = s⃗′ * s′_inv  # unit tangent vector (properly normalised) -- see also UnitTangent
    # Note: the +1 coefficient is required for energy conservation.
    # It is required so that the resulting energy follows Hamilton's equation (at least for
    # the case of a vortex ring), and it has been verified in many cases that it improves
    # the effective energy conservation.
    2 * (log(2 * δ / a) + 1 - Δ) * t̂  # note: prefactor = Γ/4π (hence the 2 in front)
end

## ========================================================================================== ##
## Computation of local integral when lia_segment_fraction < 1

function _add_local_integrals!(
        fields::NamedTuple, fs::VectorOfFilaments, cache::BiotSavartCache;
        lia_segment_fraction
    )
    lia_segment_fraction === nothing && return fields
    chunks = FilamentChunkIterator(fs)
    @sync for chunk in chunks
        isempty(chunk) && continue  # don't spawn a task if it will do no work
        Threads.@spawn for i in chunk
            prev_indices = firstindex(fs):(first(chunk) - 1)  # filament indices given to all previous chunks
            n = count_nodes(view(fs, prev_indices))  # we will start writing at index n + 1
            fields_i = map(vs -> vs[i], fields)
            n = _add_local_integrals!(fields_i, fs[i], n, cache; lia_segment_fraction)::Int
        end
    end
    fields
end

@inline function _compute_quadrature_data(f::ClosedFilament, i::Int, quad::StaticSizeQuadrature, lims::NTuple{2, Real})
    Nq = length(quad)
    ts = Filaments.knots(f)
    ζs, ws = quadrature(quad)
    a, b = lims
    δ = b - a
    Δt = @inbounds (ts[i + 1] - ts[i]) * δ
    ys = ntuple(q -> a + δ * ζs[q], Val(Nq))  # rescale and translate quadrature points to integrate within wanted limits
    s⃗_quad = ntuple(Val(Nq)) do q
        @inline
        @inbounds Tuple(f(i, ys[q]))
    end
    qs⃗′_quad = ntuple(Val(Nq)) do q
        @inline
        @inbounds Tuple(ws[q] * Δt * f(i, ys[q], Derivative(1)))
    end
    s⃗_quad, qs⃗′_quad
end

@inline function _transpose_tuples(xs::NTuple{N, NTuple{M, T}}) where {N, M, T}
    ntuple(Val(M)) do m
        @inline
        ntuple(n -> @inbounds(xs[n][m]), Val(N))::NTuple{N, T}
    end
end

function _add_local_integrals!(
        fields::NamedTuple{Names}, f::ClosedFilament, n::Int, cache::BiotSavartCache;
        lia_segment_fraction
    ) where {Names}
    lia_segment_fraction === nothing && return fields
    (; nodes) = cache.pointdata
    params = cache.params.common
    (; Γ, Ls, quad_near_singularity) = params
    T = typeof(Γ)
    lims = nonlia_integration_limits(lia_segment_fraction)
    Xs = Filaments.nodes(f)
    N = length(Ls)
    Lhs = map(L -> L / 2, Ls)
    prefactor = Γ / T(4π)
    quantities = NamedTuple{Names}(possible_output_fields())  # e.g. (velocity = Velocity(),)
    @inbounds for i in eachindex(Xs)
        n += 1
        # x⃗ = Xs[i]
        x⃗ = nodes[n]  # this should be equal to Xs[i] up to a multiple of the domain period (`nodes` is always in [0, L], unlike `Xs`)
        # Compute quadrature nodes and weights
        qdata_left = _compute_quadrature_data(f, i - 1, quad_near_singularity, lims[1])
        qdata_right = _compute_quadrature_data(f, i, quad_near_singularity, lims[2])
        s⃗_quad_unfolded = _transpose_tuples((qdata_left[1]..., qdata_right[1]...))::NTuple{N}  # one tuple per dimension
        qs⃗′_quad = _transpose_tuples((qdata_left[2]..., qdata_right[2]...))::NTuple{N}
        # Use explicit SIMD. Not sure this improves performance a lot, but it shouldn't hurt. Note
        # that the SIMD width is the total number of quadrature nodes (= 2 * length(quad_near_singularity)).
        r⃗s_simd = ntuple(Val(N)) do d
            @inline
            s_d = SIMD.Vec(map(x -> Filaments.fold_coordinates_periodic(x, Ls[d]), s⃗_quad_unfolded[d]))
            _simd_deperiodise_separation_folded(x⃗[d] - s_d, Ls[d], Lhs[d])
        end
        qs⃗′_simd = map(Vec, qs⃗′_quad)
        r²s_simd = sum(abs2, r⃗s_simd)
        rs = sqrt(r²s_simd)
        rs_inv = inv(rs)
        foreach(values(fields), values(quantities)) do vs, quantity
            @inline
            δu⃗_simd = full_integrand(quantity, rs_inv, qs⃗′_simd, r⃗s_simd)
            δu⃗_data = map(δu⃗_simd) do component
                @inline
                sum(component)  # reduction operation: sum the W elements
            end
            @inbounds vs[i] = vs[i] + prefactor * Vec3(δu⃗_data)
        end
    end
    n
end

function add_local_integrals!(
        fields::NamedTuple, cache::BiotSavartCache, fs::VectorOfFilaments,
    )
    (; lia_segment_fraction,) = cache.params.shortrange
    lia_segment_fraction === nothing && return fields  # nothing to do
    _add_local_integrals!(fields, fs, cache; lia_segment_fraction)
    fields
end
