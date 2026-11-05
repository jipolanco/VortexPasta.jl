function add_pair_interactions_2D!(outputs::NamedTuple{Names}, cache::ShortRangeCache) where {Names}
    (; splitting, pointdata, params) = cache
    (; points, signs) = pointdata
    (; Γ, Ls) = params.common
    ka_backend = KA.get_backend(cache)
    T = typeof(Γ)
    rcut² = params.rcut^2
    prefactor = Γ / T(2π)
    Lhs = map(L -> L / 2, Ls)
    quantities = NamedTuple{Names}(possible_output_fields())  # e.g. (velocity = Velocity(),)

    # We assume both source and destination points have already been folded into the main periodic cell.
    # (This is done in add_point_charges!)
    foreach_pair(cache; folded = Val(true)) do x⃗_dst, i_dst, j_src
        # @assert x⃗_dst == points[i_dst]
        if i_dst == j_src
            return nothing  # skip self-interactions
        end

        x⃗_src = @inbounds points[j_src]
        sign = @inbounds to_real(signs[j_src])
        prefactor_signed = sign * prefactor

        N = length(Ls)  # = 2
        r⃗ = ntuple(Val(N)) do d
            @inline
            local r = @inbounds x⃗_dst[d] - x⃗_src[d]
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
            a, b = weights_shortrange_nosimd(ka_backend, splitting, r²)
            args = (prefactor_signed * a, prefactor_signed * b, r², r⃗)
            map(values(outputs), values(quantities)) do vs, quantity
                @inline
                local V = eltype(vs)  # Vec2{T} or T
                @inbounds vs[i_dst] += V(short_range_contribution(quantity, args...))
                nothing
            end
        end

        nothing
    end
end

@inline short_range_contribution(::Streamfunction, a, b, r², r⃗::NTuple{2}) = a
@inline short_range_contribution(::Velocity, a, b, r², r⃗::NTuple{2}) = (b / r²) .* (-r⃗[2], +r⃗[1])
