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
    foff = end_to_end_offset(f)

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
    f2 = change_offset(f, foff - p⃗)

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
    if !iszero(foff)
        istart = lastindex(f2) - i + 1
        for k ∈ istart:lastindex(f2)
            f2[k] = f2[k] + foff
        end
    end

    f1, f2
end

"""
    merge!(f::ClosedFilament, g::ClosedFilament, i::Int, j::Int; p⃗ = Vec3(0, 0, 0)) -> ClosedFilament

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
    @debug "Merging:" f[i] g[j]
    Nf, Ng = length(f), length(g)
    is_shift = (i + 1):lastindex(f)

    foff = end_to_end_offset(f)
    goff = end_to_end_offset(g)
    offset_total = foff + goff

    resize!(f, Nf + Ng)

    # Shift second half of `f` nodes (i.e. f[is_shift]) to the end of `f`.
    l = lastindex(f) + 1
    for k ∈ reverse(is_shift)
        f[l -= 1] = f[k]
    end
    @assert l == length(f) - length(is_shift) + 1

    if !iszero(goff)
        for k ∈ firstindex(f):i
            f[k] = f[k] - goff
        end
    end

    # Copy first half of `g` nodes (i.e. g[begin:j]).
    for k ∈ j:-1:firstindex(g)
        f[l -= 1] = g[k] - p⃗
    end

    # Copy second half of `g` nodes (i.e. g[j + 1:end]).
    u⃗ = -(p⃗ + goff)
    for k ∈ lastindex(g):-1:(j + 1)
        f[l -= 1] = g[k] + u⃗
    end
    @assert l == i + 1

    change_offset(f, offset_total)
end
