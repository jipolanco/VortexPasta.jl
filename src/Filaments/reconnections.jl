using ..BasicTypes: Zero

"""
    ReconnectionCriterion

Abstract type describing a criterion for filament reconnections.

Implemented reconnection criteria include:

- [`NoReconnections`](@ref): disables reconnections;

- [`BasedOnDistance`](@ref): reconnects filament segments which are closer than a critical distance.
"""
abstract type ReconnectionCriterion end

"""
    distance(crit::ReconnectionCriterion) -> Real

Return the critical distance associated to the reconnection criterion.
"""
function distance end

"""
    NoReconnections <: ReconnectionCriterion

Used to disable filament reconnections.
"""
struct NoReconnections <: ReconnectionCriterion end

distance(::NoReconnections) = Zero()

"""
    BasedOnDistance <: ReconnectionCriterion
    BasedOnDistance(d_crit; cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

The keyword argument `cos_max` disables the reconnection of nearly parallel segments.
The default value `cos_max = 0.97` disables reconnections when the angle between lines
is ``θ < \\arccos(0.97) ≈ 14°``.
"""
struct BasedOnDistance <: ReconnectionCriterion
    dist    :: Float64
    cos_max :: Float64

    function BasedOnDistance(dist; cos_max = 0.97)
        new(dist, cos_max)
    end
end

distance(c::BasedOnDistance) = c.dist

function should_reconnect(
        c::BasedOnDistance, fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods,
    )
    (; dist_sq, cos_max_sq,) = c

    (; d⃗, ζx, ζy,) = find_min_distance(fx, fy, i, j; periods)
    d² = sum(abs2, d⃗)
    d² > dist_sq && return false  # don't reconnect

    X′ = fx(i, ζx, Derivative(1))
    Y′ = fy(j, ζy, Derivative(1))

    cos² = (X′ ⋅ Y′)^2 / (sum(abs2, X′) * sum(abs2, Y′))
    cos² < cos_max_sq
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

    inds_i = istart:lastindex(segments(f))
    periods_half = map(L -> L / 2, periods)

    # Segments to be "compared" with segment i.
    # We don't want to compare with direct neighbours, since they will often
    # pass the first filter and fail the second one (which is more expensive to compute).
    inds_j = (first(inds_i) + 2):(last(inds_i) - 1)  # avoid indices (i - 1, i, i + 1)

    for i ∈ inds_i
        isempty(inds_j) && break
        x⃗_mid = (f[i] + f[i + 1]) ./ 2
        y⃗_b = f[first(inds_j)]
        r⃗_b = deperiodise_separation(y⃗_b - x⃗_mid, periods, periods_half)
        is_outside_range_b = any(>(d_cut), r⃗_b)  # should be cheaper than computing the squared vector norm
        for j ∈ inds_j
            y⃗_a, y⃗_b = y⃗_b, f[j + 1]
            r⃗_a = deperiodise_separation(y⃗_a - x⃗_mid, periods, periods_half)
            is_outside_range_a = any(>(d_cut), r⃗_a)

            if is_outside_range_a && is_outside_range_b
                # Skip this segment if its two limits are too far from x⃗_mid.
                continue
            end

            # The current segment passed the first filter and is a candidate for reconnection.
            should_reconnect(crit, f, f, i, j; periods) || continue

            # Split filament into 2
            f₁, f₂ = split!(f, i, j)
            update_coefficients!(f₁)
            update_coefficients!(f₂)
            @assert f₁ === f
            push!(fs_new, f₂)

            # Possibly perform reconnections on each subfilament.
            # In the second case, the `istart` is to save some time by skipping
            # segment pairs which were already verified. This requires the
            # nodes in each subfilament to be sorted in a specific manner, and
            # may fail if the split! function is modified.
            reconnect_self!(crit, f₁, fs_new; periods)
            reconnect_self!(crit, f₂, fs_new; periods, istart = i + 1)
        end
        inds_j = (i + 3):last(inds_i)  # for next iteration
    end

    fs_new
end

"""
    split!(f::ClosedFilament, i::Int, j::Int) -> (f₁, f₂)

Split closed filament into two filaments.

Assuming `j > i`, the resulting filaments are respectively composed of nodes
`f[i + 1:j]` and `f[(begin:i) ∪ (j + 1:end)]`.

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
    for k ∈ firstindex(f):i
        f2[l += 1] = f[k]
    end
    for k ∈ (j + 1):lastindex(f)
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
