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
  Note that the angle ``θ`` is signed (it takes values in ``[-1, 1]``).
  Negative angles mean that the segments are antiparallel, and in this case reconnections are always performed.
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
        crit::ReconnectionCriterion,
        f::AbstractFilament,
        flist::AbstractVector{<:AbstractFilament};
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) -> Union{Nothing, AbstractFilament}

Attempt to reconnect filament `f` with itself.

Note that reconnections of a filament with itself produce new filaments.
Newly created filaments will be appended to the `flist` container.

This function returns `nothing` if no reconnection happened.

Otherwise, if a reconnection happened, this function returns one of the
resulting filaments, `f₁`. The other (one or more) resulting filaments are
appended to `flist`.

For example, if the filament `f` self-reconnects onto 4 filaments, this
function returns one of these filaments, and the other 3 are appended to
`flist`. This is useful if one wants to replace the original filament `f` by
one of the resulting filaments.

In a periodic domain, one should also pass the optional `periods` argument.
"""
function reconnect_self!(
        crit::ReconnectionCriterion, f::F,
        fs_new::AbstractVector{F} = F[];
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
        istart = firstindex(segments(f)),
    ) where {F <: AbstractFilament}
    d_crit = distance(crit)
    d_crit === Zero() && return nothing  # reconnections are disabled

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
            f₁, f₂ = split!(f, i, j; p⃗ = info.p⃗)

            # Update coefficients and possibly perform reconnections on each subfilament.
            # In the first case, the `istart` is to save some time by skipping
            # segment pairs which were already verified. This requires the
            # nodes in each subfilament to be sorted in a specific manner, and
            # may fail if the split! function is modified.
            if check_nodes(Bool, f₂)  # skip if coefficients can't be computed, typically if the number of nodes is too small (< 3 for cubic splines)
                update_coefficients!(f₂)
                g₂ = reconnect_self!(crit, f₂, fs_new; periods, istart = i + 1)
                push!(fs_new, something(g₂, f₂))  # push f₂, or its replacement if f₂ itself reconnected
            end

            if check_nodes(Bool, f₁)
                update_coefficients!(f₁)
                g₁ = reconnect_self!(crit, f₁, fs_new; periods)
                return something(g₁, f₁)  # return f₁, or its replacement if f₁ itself reconnected
            end

            return f₁  # we can stop iterating here, since the filament `f` doesn't exist anymore
        end
        inds_j = (i + 3):last(inds_i)  # for next iteration
    end

    nothing
end

"""
    reconnect_other!(
        crit::ReconnectionCriterion, f::AbstractFilament, g::AbstractFilament;
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    ) -> Union{Nothing, AbstractFilament}

Attempt to reconnect filaments `f` and `g`.

The two filaments cannot be the same. To reconnect a filament with itself, see [`reconnect_self!`](@ref).

This function allows at most a single reconnection between the two filaments.
If a reconnection happens, the two filaments merge into one, and the resulting filament is returned.
The original filaments `f` and `g` can be discarded (in particular, `f` is modified internally).

Returns the merged filament a reconnection happened, `nothing` otherwise.
"""
function reconnect_other!(
        crit::ReconnectionCriterion, f::AbstractFilament, g::AbstractFilament;
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
    )
    @assert f !== g

    d_crit = distance(crit)
    d_crit === Zero() && return nothing  # reconnections are disabled

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

            h = merge!(f, g, i, j; p⃗ = info.p⃗)  # filaments are merged onto `h`
            update_coefficients!(h)
            return h
        end
    end

    nothing
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

This function calls [`reconnect_other!`](@ref) on all filament pairs, and then
[`reconnect_self!`](@ref) on all filaments.
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
            h = reconnect_other!(crit, f, g; periods)
            h === nothing && continue

            # The two filaments were merged into `h`, and filaments `f` and `g` can be removed.
            fs[i] = h
            callback(h, i, :modified)
            popat!(fs, j)
            callback(g, j, :removed)
            j -= 1
            jlast -= 1
            ilast -= 1
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
        f₁ = reconnect_self!(crit, f, fs; periods)
        f₁ === nothing && continue  # there were no reconnections

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

        # Now replace the old filament `f` by filament `f₁` (which is one of
        # the filaments resulting from the reconnection).
        if check_nodes(Bool, f₁)
            fs[i] = f₁
            callback(f₁, i, :modified)
        else
            popat!(fs, i)
            callback(f₁, i, :removed)
            i -= 1
            ilast -= 1
        end
    end

    fs
end

reconnect!(crit::ReconnectionCriterion, args...; kws...) =
    reconnect!(Returns(nothing), crit, args...; kws...)

"""
    split!(f::ClosedFilament, i::Int, j::Int; p⃗ = Vec3(0, 0, 0)) -> (f₁, f₂)

Split closed filament into two filaments.

Assuming `j > i`, the resulting filaments are respectively composed of nodes
`f[i + 1:j]` and `f[(j + 1:end) ∪ (begin:i)]`.

In practice, a split makes sense when the nodes `f[i]` and `f[j] - p⃗` are spatially "close".
Here `p⃗` is an optional offset which usually takes into account domain periodicity
(see also [`find_min_distance`](@ref)).

Note that this function modifies the filament `f`, which should then be discarded.

One should generally call [`update_coefficients!`](@ref) on both filaments after a split.
"""
function split!(f::ClosedFilament, i::Int, j::Int; p⃗ = Vec3(0, 0, 0))
    i ≤ j || return split!(f, j, i)
    @debug "Splitting:" f[i] f[j] p⃗

    n1 = j - i
    n2 = length(f) - n1

    # Copy values onto f1
    f1 = change_offset(similar(f, n1), p⃗)
    l = firstindex(f1) - 1
    for k ∈ (i + 1):j
        f1[l += 1] = f[k]
    end
    @assert l == n1

    # Try to reuse `f` to reduce allocations.
    # For some reason here the offset must be relative to the original one, which was not
    # the case for f1.
    f2 = change_offset(f, f.Xoffset - p⃗)

    # Fill f2.
    # We want nodes to be in the order f[(j + 1:end) ∪ (begin:i)], which is more difficult
    # to do in-place compared to the order f[(begin:i) ∪ (j + 1:end)].
    # Only the first order is valid when there's a periodic offset p⃗ (since the offset is
    # between f[i] and f[j + 1], and so these should be at the limits of the new filament).

    # We start by doing a circular shift of the nodes of `f` (which are also the nodes of
    # `f2`), so that f[j + 1] is moved to the beginning.
    # We use the trick described in https://stackoverflow.com/a/44658599, which requires
    # three in-place reversals.
    Xs = nodes(f)
    @assert Xs === nodes(f2)
    shift = j + 1 - firstindex(Xs)
    reverse!(Xs)
    reverse!(Xs, firstindex(Xs), lastindex(Xs) - shift)
    reverse!(Xs, lastindex(Xs) - shift + 1, lastindex(Xs))

    # Now, the nodes of f2 are simply the `n2` first elements of Xs, so we just
    # need to resize f2 and we're (almost) done.
    resize!(f2, n2)

    # Finally, we need to take into account the offset of the original filament `f`, to make
    # sure that we don't have a jump between the ranges (j + 1:end) and (begin:i).
    off = f.Xoffset
    if !iszero(off)
        istart = lastindex(f2) - i + 1
        for k ∈ istart:lastindex(f2)
            f2[k] = f2[k] + off
        end
    end

    f1, f2
end

"""
    merge!(f::ClosedFilament, g::ClosedFilament, i::Int, j::Int; p⃗ = Vec3(0, 0, 0))

Merge two closed filaments into one.

The resulting filament is composed of nodes:

    f[begin:i] ∪ {g[(j + 1:end) ∪ (begin:j)] - p⃗} ∪ f[i + 1:end]

Here `p⃗` is an optional offset which usually takes into account domain periodicity
(see also [`find_min_distance`](@ref)).

This function returns a new filament `h` which may share memory with `f`.
The filament `g` is not modified.

One should generally call [`update_coefficients!`](@ref) on the returned filament after merging.
"""
function merge!(f::ClosedFilament, g::ClosedFilament, i::Int, j::Int; p⃗::Vec3 = zero(eltype(f)))
    @debug "Merging:" f[i] g[j] p⃗ f.Xoffset g.Xoffset
    Nf, Ng = length(f), length(g)
    is_shift = (i + 1):lastindex(f)

    resize!(f, Nf + Ng)
    foff = f.Xoffset
    goff = g.Xoffset

    # Shift second half of `f` nodes (i.e. f[is_shift]) to the end of `f`.
    l = lastindex(f) + 1
    for k ∈ reverse(is_shift)
        f[l -= 1] = f[k]
    end
    @assert l == length(f) - length(is_shift) + 1
    for k ∈ firstindex(f):i
        f[k] = f[k] + foff
    end

    # Copy first half of `g` nodes (i.e. g[begin:j]).
    for k ∈ j:-1:firstindex(g)
        f[l -= 1] = g[k] - p⃗
    end

    # Copy second half of `g` nodes (i.e. g[j + 1:end]).
    u⃗ = p⃗ + goff
    for k ∈ lastindex(g):-1:(j + 1)
        f[l -= 1] = g[k] - u⃗
    end
    @assert l == i + 1

    off = foff + goff
    change_offset(f, off)
end
