using ..BasicTypes: Zero, Infinity
using LinearAlgebra: ⋅, norm

"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`BasedOnDistance`](@ref): reconnects filament segments which are closer than a critical distance.
"""
abstract type ReconnectionCriterion end

"""
    should_reconnect(
        c::ReconnectionCriterion,
        fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
    ) -> Union{Nothing, NamedTuple}

Check whether two filaments should reconnect according to the chosen criterion.

Checks for a possible reconnection between filament segments `fx[i:i+1]` and `fy[j:j+1]`.

If the filament segments should reconnect, this function returns a `NamedTuple`
with reconnection information, which includes in particular all the fields
returned by [`find_min_distance`](@ref).

Otherwise, returns `nothing` if the filament segments should not reconnect.
"""
function should_reconnect end

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = Zero()
should_reconnect(::NoReconnections, args...; kws...) = nothing

"""
    BasedOnDistance <: ReconnectionCriterion
    BasedOnDistance(d_crit; decrease_length = true, cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

# Optional keyword arguments

- `decrease_length`: if `true` (default), a reconnection will only be performed
  if it will decrease the total filament length. Since, for vortices, the total
  energy is roughly related to the vortex length, this means that reconnections
  should always tend to dissipate energy.

- `cos_max`: allows to disable reconnections of nearly parallel segments. Two segments
  are considered to be "nearly parallel" if `cos(θ) > cos_max`.
  The default value `cos_max = 0.97` disables reconnections when the angle between lines
  is ``θ < \\arccos(0.97) ≈ 14°``.
"""
struct BasedOnDistance <: ReconnectionCriterion
    dist       :: Float64
    dist_sq    :: Float64
    cos_max    :: Float64
    cos_max_sq :: Float64
    decrease_length :: Bool

    function BasedOnDistance(dist; cos_max = 0.97, decrease_length = true)
        new(dist, dist^2, cos_max, cos_max^2, decrease_length)
    end
end

distance(c::BasedOnDistance) = c.dist

function should_reconnect(
        c::BasedOnDistance, fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
    )
    (; dist_sq, cos_max_sq, decrease_length,) = c

    min_dist = find_min_distance(fx, fy, i, j; periods)
    (; d⃗, p⃗, ζx, ζy,) = min_dist
    d² = sum(abs2, d⃗)
    d² > dist_sq && return nothing  # don't reconnect

    # Make sure that reconnections reduce the total length (makes sense energetically for vortices).
    if decrease_length
        length_before = norm(fx[i + 1] - fx[i]) + norm(fy[j + 1] - fy[j])
        length_after = norm(fy[j + 1] - fx[i] - p⃗) + norm(fx[i + 1] - fy[j] + p⃗)
        length_after > length_before && return nothing
    end

    X′ = fx(i, ζx, Derivative(1))
    Y′ = fy(j, ζy, Derivative(1))

    success = min_dist  # for now, only return the output of find_min_distance if segments should reconnect

    xy = X′ ⋅ Y′
    xy < 0 && return success  # always reconnect antiparallel vortices

    cos² = (xy * xy) / (sum(abs2, X′) * sum(abs2, Y′))
    cos² < cos_max_sq ? success : nothing
end

"""
    reconnect_self!(
        crit::ReconnectionCriterion, f::AbstractFilament,
        [fs::AbstractVector{<:AbstractFilament}];
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) -> fs_added

Attempt to reconnect filament `f` with itself.

Note that reconnections of a filament with itself produce new filaments.
By default, this function allocates a new vector which will contain the newly created filaments.
Optionally, one may pass an existent filament container `fs` to which new filaments will be appended.

In both cases, this function returns a vector view `fs_added` containing the newly created filaments.
This means that one can get the number of added filaments by just doing `length(fs_added)`.

In a periodic domain, one should also pass the optional `periods` argument.
"""
function reconnect_self!(
        crit::ReconnectionCriterion, f::F,
        fs_new::AbstractVector{F} = F[];
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
        istart = firstindex(segments(f)),
    ) where {F <: AbstractFilament}
    d_crit = distance(crit)
    d_crit === Zero() && return fs_new  # reconnections are disabled

    # This cutoff distance serves as a first (coarse) filter.
    # It is larger than the critical distance to take into account the fact
    # that the coarse filter doesn't use the *minimal* distance between
    # segments, but only a rough approximation.
    d_cut = 2 * d_crit
    d_cut_sq = d_cut^2

    inds_i = istart:lastindex(segments(f))
    periods_half = map(L -> L / 2, periods)

    # Segments to be "compared" with segment i.
    # We don't want to compare with direct neighbours, since they will often
    # pass the first filter and fail the second one (which is more expensive to compute).
    inds_j = (first(inds_i) + 2):(last(inds_i) - 1)  # avoid indices (i - 1, i, i + 1)

    for i ∈ inds_i
        isempty(inds_j) && break
        x⃗_mid = (f[i] + f[i + 1]) ./ 2
        r⃗_b = let
            y⃗ = f[first(inds_j)]
            deperiodise_separation(y⃗ - x⃗_mid, periods, periods_half)
        end
        is_outside_range_b = any(>(d_cut), r⃗_b)  # should be cheaper than computing the squared vector norm
        for j ∈ inds_j
            r⃗_a = r⃗_b
            r⃗_b = let
                y⃗ = f[j + 1]
                deperiodise_separation(y⃗ - x⃗_mid, periods, periods_half)
            end
            is_outside_range_a = is_outside_range_b
            is_outside_range_b = any(>(d_cut), r⃗_b)

            if is_outside_range_a && is_outside_range_b
                # Skip this segment if its two limits are too far from x⃗_mid.
                continue
            end

            # Second (slightly finer) filter: look at the actual distances.
            if sum(abs2, r⃗_a) > d_cut_sq || sum(abs2, r⃗_b) > d_cut_sq
                continue
            end

            # The current segment passed the first two filters and is a candidate for reconnection.
            info = should_reconnect(crit, f, f, i, j; periods)
            info === nothing && continue

            # Split filament into 2
            f₁, f₂ = split!(f, i, j)
            @assert f₁ === f
            push!(fs_new, f₂)

            # Update coefficients and possibly perform reconnections on each subfilament.
            # In the second case, the `istart` is to save some time by skipping
            # segment pairs which were already verified. This requires the
            # nodes in each subfilament to be sorted in a specific manner, and
            # may fail if the split! function is modified.
            if check_nodes(Bool, f₁)  # skip if coefficients can't be computed, typically if the number of nodes is too small (< 3 for cubic splines)
                update_coefficients!(f₁)
                reconnect_self!(crit, f₁, fs_new; periods)
            end

            if check_nodes(Bool, f₂)
                update_coefficients!(f₂)
                reconnect_self!(crit, f₂, fs_new; periods, istart = i + 1)
            end

            return fs_new  # I don't need to continue iterating on this filament (which is the same as f₁)
        end
        inds_j = (i + 3):last(inds_i)  # for next iteration
    end

    fs_new
end

"""
    reconnect_other!(
        crit::ReconnectionCriterion, f::AbstractFilament, g::AbstractFilament;
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )

Attempt to reconnect filaments `f` and `g`.

The two filaments cannot be the same. To reconnect a filament with itself, see [`reconnect_self!`](@ref).

This function allows at most a single reconnection between the two vortices.
If a reconnection happens, the two vortices merge into one. Vortex `f` is
updated with the result, while vortex `g` can be discarded.

Returns `true` if a reconnection happened, `false` otherwise.
"""
function reconnect_other!(
        crit::ReconnectionCriterion, f::AbstractFilament, g::AbstractFilament;
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )
    @assert f !== g

    d_crit = distance(crit)
    d_crit === Zero() && return false  # reconnections are disabled

    # The following is very similar to `reconnect_self!`
    d_cut = 2 * d_crit
    d_cut_sq = d_cut^2
    periods_half = map(L -> L / 2, periods)

    inds_i = eachindex(segments(f))
    inds_j = eachindex(segments(g))

    for i ∈ inds_i
        x⃗_mid = (f[i] + f[i + 1]) ./ 2
        r⃗_b = let
            y⃗ = g[first(inds_j)]
            deperiodise_separation(y⃗ - x⃗_mid, periods, periods_half)
        end
        is_outside_range_b = any(>(d_cut), r⃗_b)  # should be cheaper than computing the squared vector norm
        for j ∈ inds_j
            r⃗_a = r⃗_b
            r⃗_b = let
                y⃗ = g[j + 1]
                deperiodise_separation(y⃗ - x⃗_mid, periods, periods_half)
            end
            is_outside_range_a = is_outside_range_b
            is_outside_range_b = any(>(d_cut), r⃗_b)

            if is_outside_range_a && is_outside_range_b
                # Skip this segment if its two limits are too far from x⃗_mid.
                continue
            end

            # Second (slightly finer) filter: look at the actual distances.
            if sum(abs2, r⃗_a) > d_cut_sq || sum(abs2, r⃗_b) > d_cut_sq
                continue
            end

            # The current segment passed the first two filters and is a candidate for reconnection.
            info =  should_reconnect(crit, f, g, i, j; periods)
            info === nothing && continue

            merge!(f, g, i, j; p⃗ = info.p⃗)  # filaments are merged onto `f`
            update_coefficients!(f)
            return true
        end
    end

    false
end

"""
    reconnect!(
        [callback::Function],
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament};
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )

Perform filament reconnections according to chosen criterion.

Note that, when a filament self-reconnects, this creates new filaments, which
are appended at the end of `fs`.

Moreover, this function will remove reconnected filaments if their number of nodes is too small
(typically ``< 3``, see [`check_nodes`](@ref)).

## Callback function

Optionally, one may pass a callback which will be called whenever the vector of
filaments `fs` is modified. Its signature must be the following:

    callback(f::AbstractFilament, i::Int, mode::Symbol)

where `f` is the modified filament, `i` is its index in `fs`, and `mode` is one of:

- `:modified` if the filament `fs[i]` was modified;
- `:appended` if the filament was appended at the end of `fs` (at index `i`);
- `:removed` if the filament previously located at index `i` was removed.

See also [`reconnect_self!`](@ref).
"""
function reconnect!(
        callback::F,
        crit::ReconnectionCriterion,
        fs::AbstractVector{<:AbstractFilament};
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) where {F <: Function}

    # 1. Reconnect filaments with each other.
    i = firstindex(fs) - 1
    ilast = lastindex(fs)
    while i < ilast
        i += 1
        f = fs[i]
        j = i
        jlast = lastindex(fs)
        while j < jlast
            j += 1
            g = fs[j]
            if reconnect_other!(crit, f, g; periods)
                # The two filaments were merged into `f`, and the filament `g` can be removed.
                callback(f, i, :modified)
                popat!(fs, j)
                callback(g, j, :removed)
                j -= 1
                jlast -= 1
                ilast -= 1
            end
        end
    end

    # 2. Reconnect filaments with themselves.
    # This needs to be done after reconnecting with each other.
    # In the specific case of two vortices reconnecting at two separate locations, this ensures that:
    #  (i)  In step 1, the two vortices reconnect at one of the locations forming one vortex in step 1.
    #  (ii) In step 2, the "big" vortex generated in step 1 self-reconnects, ending up with two vortices.
    i = firstindex(fs) - 1
    ilast = lastindex(fs)

    while i < ilast  # make sure we don't include appended filaments in the iteration
        i += 1
        f = fs[i]
        n_old = lastindex(fs)
        reconnect_self!(crit, f, fs; periods)

        if lastindex(fs) === n_old  # there were no reconnections
            continue
        end

        # First check appended filaments, and remove them if they don't have enough nodes (typically < 3).
        j = n_old
        while j < lastindex(fs)
            j += 1
            fj = fs[j]  # this is an appended filament
            if check_nodes(Bool, fj)
                callback(fj, j, :appended)
            else
                popat!(fs, j)  # remove the filament
                j -= 1
            end
        end

        # Now check the modified filament f = fs[i].
        if check_nodes(Bool, f)
            callback(f, i, :modified)
        else
            popat!(fs, i)
            callback(f, i, :removed)
            i -= 1
            ilast -= 1
        end
    end

    fs
end

reconnect!(crit::ReconnectionCriterion, args...; kws...) =
    reconnect!(Returns(nothing), crit, args...; kws...)

"""
    split!(f::ClosedFilament, i::Int, j::Int) -> (f₁, f₂)

Split closed filament into two filaments.

Assuming `j > i`, the resulting filaments are respectively composed of nodes
`f[i + 1:j]` and `f[(j + 1:end) ∪ (begin:i)]`.

In practice, a split makes sense when the nodes `f[i]` and `f[j]` are spatially "close".

Note that the filament `f` is modified by this function, and is returned as the filament `f₁`.

One should generally call [`update_coefficients!`](@ref) on both filaments after a split.
"""
function split!(f::ClosedFilament, i::Int, j::Int)
    i ≤ j || return split!(f, j, i)

    n1 = j - i
    n2 = length(f) - n1

    f1 = f
    f2 = similar(f, n2)

    # Fill f2
    l = firstindex(f2) - 1
    for k ∈ (j + 1):lastindex(f)
        f2[l += 1] = f[k]
    end
    for k ∈ firstindex(f):i
        f2[l += 1] = f[k]
    end
    @assert l == n2

    # Shift values and then resize f1
    l = firstindex(f1) - 1
    for k ∈ (i + 1):j
        f1[l += 1] = f[k]
    end
    @assert l == n1
    resize!(f1, n1)

    f1, f2
end

"""
    merge!(f::ClosedFilament, g::ClosedFilament, i::Int, j::Int; p⃗ = Vec3(0, 0, 0))

Merge two closed filaments into one.

The resulting filament is composed of nodes:

    f[begin:i] ∪ {g[(j + 1:end) ∪ (begin:j)] - p⃗} ∪ f[i + 1:end]

Here `p⃗` is an optional offset which usually takes into account domain periodicity
(see also [`find_min_distance`](@ref)).

This function modifies (and returns) the filament `f`, which contains the
merged filament at output. The filament `g` is not modified.

One should generally call [`update_coefficients!`](@ref) on `f` after merging.
"""
function merge!(f::ClosedFilament, g::ClosedFilament, i::Int, j::Int; p⃗::Vec3 = zero(eltype(f)))
    Nf, Ng = length(f), length(g)
    is_shift = (i + 1):lastindex(f)

    resize!(f, Nf + Ng)

    # Shift second half of `f` nodes (i.e. f[is_shift]) to the end of `f`.
    l = lastindex(f) + 1
    for k ∈ reverse(is_shift)
        f[l -= 1] = f[k]
    end
    @assert l == length(f) - length(is_shift) + 1

    # Copy first half of `g` nodes (i.e. g[begin:j]).
    for k ∈ j:-1:firstindex(g)
        f[l -= 1] = g[k] - p⃗
    end

    # Copy second half of `g` nodes (i.e. g[j + 1:end]).
    for k ∈ lastindex(g):-1:(j + 1)
        f[l -= 1] = g[k] - p⃗
    end
    @assert l == i + 1

    f
end
