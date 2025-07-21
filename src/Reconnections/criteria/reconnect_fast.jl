using ..CellLists: CellLists, PeriodicCellList
using LinearAlgebra: norm, ⋅

"""
    ReconnectFast <: ReconnectionCriterion
    ReconnectFast(d_crit; decrease_length = true, cos_max = 0.97)

Reconnects filament segments which are at a distance `d < d_crit`.

This criterion is similar to [`ReconnectBasedOnDistance`](@ref) but with important differences:

- filament segments are considered as straight lines (which may be less accurate);

- all reconnections are made at once, i.e. each filament will be reconnected multiple times if needed.
  In other words, a single pass should be needed to reconnect all segments satisfying the chosen criterion.

This criterion only supports periodic domains.

See [`ReconnectBasedOnDistance`](@ref) for details on criterion parameters.
"""
struct ReconnectFast <: ReconnectionCriterion
    dist       :: Float64
    dist_sq    :: Float64
    cos_max    :: Float64
    cos_max_sq :: Float64
    use_velocity    :: Bool  # this is currently ignored (TODO: implement?)
    decrease_length :: Bool
    function ReconnectFast(dist; cos_max = 0.97, use_velocity = false, decrease_length = true)
        new(dist, dist^2, cos_max, cos_max^2, use_velocity, decrease_length)
    end
end

distance(c::ReconnectFast) = c.dist
max_passes(c::ReconnectFast) = 1

function Base.show(io::IO, c::ReconnectFast)
    print(io, "ReconnectFast($(c.dist); cos_max = $(c.cos_max), use_velocity = $(c.use_velocity), decrease_length = $(c.decrease_length))")
end

struct ReconnectFastCache{
        T <: AbstractFloat,
        Criterion <: ReconnectFast,
        Finder <: PeriodicCellList,
        Periods <: NTuple{3, T},
        PointsVector <: AbstractVector{Vec3{T}},
        IndexVector <: AbstractVector{<:Integer},
    } <: AbstractReconnectionCache
    crit::Criterion
    cl::Finder
    Ls::Periods
    nodes::PointsVector          # length Np = number of nodes (filament points)
    filament_modified::Vector{Bool}  # length Nf; whether each filament has been modified
    filament_idx::IndexVector    # length Np; filament index each node currently belongs to
    node_next::IndexVector       # length Np; values in 1:Np
    node_prev::IndexVector       # length Np; values in 1:Np
    node_reconnected::BitVector  # length Np
end

criterion(c::ReconnectFastCache) = c.crit

function _init_cache(crit::ReconnectFast, fs::AbstractVector{<:AbstractFilament}, Ls_in::Tuple)
    T = number_type(fs)
    has_nonperiodic_directions = any(L -> L === Infinity(), Ls_in)
    has_nonperiodic_directions && throw(ArgumentError("ReconnectFast only allows fully periodic domains"))
    Ls = map(T, Ls_in)    # convert number type if needed
    nodes = similar(Filaments.nodes(first(fs)), 0)
    node_next = Vector{Int}(undef, 0)
    node_prev = similar(node_next)
    filament_modified = Vector{Bool}(undef, 0)
    node_reconnected = falses(0)
    filament_idx = similar(node_next)
    d_reconnect = distance(crit)
    r_cut = T(d_reconnect)
    Lmin = T(min(Ls...))
    M = 2  # number of subdivisions
    r_cut_max = CellLists.max_cutoff_distance(M, Lmin)
    r_cut > r_cut_max && error(
        lazy"""reconnection distance is too large compared to the domain size: d_reconnect / L_min = $(d_reconnect / Lmin).
        Try a smaller reconnection distance or a larger domain."""
    )
    rs = map(_ -> r_cut, Ls)
    nsubdiv = Val(M)
    I = eltype(node_next)  # integer type
    cl = PeriodicCellList(I, rs, Ls, nsubdiv)  # a single element is a node index (in 1:Np)
    ReconnectFastCache(crit, cl, Ls, nodes, filament_modified, filament_idx, node_next, node_prev, node_reconnected)
end

function _update_cache!(cache::ReconnectFastCache, fs::AbstractVector{<:ClosedFilament}; to)
    (; cl, nodes, node_next, node_prev, filament_idx, filament_modified, node_reconnected,) = cache
    Np = sum(length, fs)  # total number of filament nodes
    Nf = length(fs)       # number of filaments
    resize!(nodes, Np)
    resize!(node_next, Np)
    resize!(node_prev, Np)
    resize!(filament_idx, Np)       # TODO: do we use this?
    resize!(filament_modified, Nf)  # TODO: do we use this?
    empty!(cl)  # remove previous elements
    sizehint!(cl, Np)
    resize!(node_reconnected, Np)
    fill!(node_reconnected, false)  # TODO: do we use this?
    n = 0
    @inbounds for (i, f) in pairs(fs)
        filament_modified[i] = false
        nstart = n + 1  # this filament currently starts at node n + 1
        xs = Filaments.nodes(f)
        for x in xs
            n += 1
            CellLists.add_element!(cl, n, x)
            filament_idx[n] = i
            nodes[n] = x
            node_prev[n] = n - 1
            node_next[n] = n + 1
        end
        # "Close" the filament
        node_prev[nstart] = n
        node_next[n] = nstart
    end
    CellLists.finalise_cells!(cl)
    cache
end

function should_reconnect(crit::ReconnectFast, nodes, i, j; Ls, node_prev, node_next)
    (; dist_sq, cos_max, decrease_length,) = crit

    # Check that we're not trying to reconnect the same or neighbouring points.
    i == j && return nothing
    @inbounds let i_next = i, j_next = j
        for _ in 1:2
            i_next = node_next[i_next]
            j == i_next && return nothing
            j_next = node_next[j_next]
            i == j_next && return nothing
        end
    end

    Lhs = map(L -> L / 2, Ls)
    x⃗ = @inbounds nodes[i]
    y⃗ = @inbounds nodes[j]
    d⃗ = deperiodise_separation(x⃗ - y⃗, Ls, Lhs)
    d² = sum(abs2, d⃗)

    if d² > dist_sq
        return nothing
    end

    # Determine best segment to reconnect.
    # For this, we want to find the nearest pair of segments containing points (i, j).
    # There are 2×2 = 4 possible combinations.
    #
    # We use deperiodise_separation just in case we're at the start or the end of an
    # infinite filament with non-zero end-to-end offset (i.e. an unclosed filament).
    # We could avoid this if we included points `begin - 1` and `end + 1` in the `nodes` vector
    # (but this might introduce other complications).
    i_prev, i_next = node_prev[i], node_next[i]
    j_prev, j_next = node_prev[j], node_next[j]
    x⃗_prev = x⃗ + deperiodise_separation(nodes[i_prev] - x⃗, Ls, Lhs)  # this is usually just nodes[i_prev]
    x⃗_next = x⃗ + deperiodise_separation(nodes[i_next] - x⃗, Ls, Lhs)
    y⃗_prev = y⃗ + deperiodise_separation(nodes[j_prev] - y⃗, Ls, Lhs)
    y⃗_next = y⃗ + deperiodise_separation(nodes[j_next] - y⃗, Ls, Lhs)

    # Segment midpoints
    x⃗_mid_prev = (x⃗_prev + x⃗) ./ 2
    x⃗_mid_next = (x⃗ + x⃗_next) ./ 2
    y⃗_mid_prev = (y⃗_prev + y⃗) ./ 2
    y⃗_mid_next = (y⃗ + y⃗_next) ./ 2

    # Four possible segment combinations
    cases = (
        (
            (i_prev, i, x⃗_prev, x⃗, x⃗_mid_prev),
            (j_prev, j, y⃗_prev, y⃗, y⃗_mid_prev),
        ),
        (
            (i_prev, i, x⃗_prev, x⃗, x⃗_mid_prev),
            (j, j_next, y⃗, y⃗_next, y⃗_mid_next),
        ),
        (
            (i, i_next, x⃗, x⃗_next, x⃗_mid_next),
            (j_prev, j, y⃗_prev, y⃗, y⃗_mid_prev),
        ),
        (
            (i, i_next, x⃗, x⃗_next, x⃗_mid_next),
            (j, j_next, y⃗, y⃗_next, y⃗_mid_next),
        ),
    )

    d²_seg_min = oftype(d², Inf)  # minimum distance between segment midpoints
    ncase_best = 0

    for (ncase, case) in pairs(cases)
        px, py = case
        x⃗_mid = px[5]
        y⃗_mid = py[5]
        d⃗_seg = deperiodise_separation(x⃗_mid - y⃗_mid, Ls, Lhs)
        d²_seg = sum(abs2, d⃗_seg)
        if d²_seg < d²_seg_min
            ncase_best = ncase
            d²_seg_min = d²_seg
        end
    end

    case_best = @inbounds cases[ncase_best]
    i⁻, i⁺, x⃗⁻, x⃗⁺, _ = case_best[1]
    j⁻, j⁺, y⃗⁻, y⃗⁺, _ = case_best[2]

    # Separation of segments after reconnection.
    # Often, of these two is exactly equal to the d² computed above.
    δa² = sum(abs2, deperiodise_separation(y⃗⁺ - x⃗⁻, Ls, Lhs))
    δb² = sum(abs2, deperiodise_separation(x⃗⁺ - y⃗⁻, Ls, Lhs))
    if min(δa², δb²) < d²
        # If one of the distances is smaller than d, it means that the original (i, j)
        # pair is not the minimum distance between nodes, and there might be a better
        # reconnection candidate within neighbouring segments. So we stop here and keep
        # looking for other candidates.
        return nothing
    end

    δa = sqrt(δa²)
    δb = sqrt(δb²)

    # Separation of segments to be reconnected (with the right orientation)
    δx⃗ = x⃗⁺ - x⃗⁻  # no need to deperiodise this one, it's already done
    δy⃗ = y⃗⁺ - y⃗⁻
    δx = norm(δx⃗)
    δy = norm(δy⃗)

    length_before = δx + δy
    length_after = δa + δb
    if decrease_length && length_after > length_before
        return nothing
    end

    info = (;
        length_before, length_after,
        is = (i⁻, i⁺),
        js = (j⁻, j⁺),
    )

    xy = δx⃗ ⋅ δy⃗
    xy < 0 && return info  # always reconnect antiparallel vortices

    cos_xy = xy / (δx * δy)

    if cos_xy < cos_max
        info
    else
        nothing
    end
end

function _reconnect_from_cache!(cache::ReconnectFastCache; to)
    (; crit, cl, nodes, node_prev, node_next, node_reconnected, filament_idx, filament_modified, Ls,) = cache
    # TODO: parallelise?
    reconnection_count = 0
    reconnection_length_loss = zero(number_type(nodes))
    @inbounds for i in eachindex(nodes)
        x⃗ = nodes[i]
        for j in CellLists.nearby_elements(cl, x⃗)
            # i == j && continue
            info = should_reconnect(crit, nodes, i, j; Ls, node_prev, node_next)
            info === nothing && continue
            # fi = filament_idx[i]  # this is the filament index before any reconnections
            # fj = filament_idx[j]
            # same_filament = fi == fj
            (; is, js, length_before, length_after,) = info

            # TODO: remove this?
            # already_reconnected = any(n -> node_reconnected[n], (is..., js...))
            # already_reconnected && continue
            # for n in (is..., js...)
            #     # @assert node_reconnected[n] == false
            #     node_reconnected[n] = true
            # end

            i⁻, i⁺ = is
            j⁻, j⁺ = js
            reconnection_count += 1
            reconnection_length_loss += length_before - length_after
            # Reconnect i⁻ -> j⁺ and j⁻ -> i⁺
            node_next[i⁻] = j⁺
            node_prev[j⁺] = i⁻
            node_next[j⁻] = i⁺
            node_prev[i⁺] = j⁻
            # if same_filament
            #     # filament_modified[fi] = true  # TODO: is this used?
            #     Nf_after += 1  # self-reconnection creates a new filament
            # else
            #     # filament_modified[fi] = true
            #     # filament_modified[fj] = true
            #     Nf_after -= 1  # two filaments were merged
            # end
        end
    end
    # @show length(filament_modified) sum(filament_modified) Nf_after
    (; reconnection_count, reconnection_length_loss,)
end

@inline function periodic_offset(x⃗::Vec3, y⃗::Vec3, Ls::Tuple, Lhs::Tuple)
    Vec3(map(periodic_offset, Tuple(x⃗), Tuple(y⃗), Ls, Lhs))
end

@inline function periodic_offset(x::Real, y::Real, L, Lh)
    pL = zero(L)
    d = y - x
    while d > Lh
        d -= L
        pL += L
    end
    while d <= -Lh
        d += L
        pL -= L
    end
    pL
end

function _reconstruct_filaments!(callback::F, fs, cache::ReconnectFastCache) where {F}
    (; nodes, node_prev, node_next, filament_idx, Ls,) = cache
    filaments_removed_count = 0
    filaments_removed_length = zero(number_type(nodes))
    Lhs = map(L -> L / 2, Ls)
    # For now, simply delete all previous filaments and recreate new ones.
    fbase = first(fs)  # this is just to make it easier to create new filaments of the right type
    method = Filaments.discretisation_method(fbase)
    for i in reverse(eachindex(fs))
        f = pop!(fs)
        callback(f, i, :removed)
    end
    @assert isempty(fs)
    reconstructed = falses(length(nodes))
    i_start = firstindex(nodes)
    while i_start <= lastindex(nodes)
        # Reconstruct a single filament
        reconstructed[i_start] = true
        x_first = nodes[i_start]
        x_prev = x_first
        xs = [x_first]
        i = node_next[i_start]
        while i != i_start
            reconstructed[i] = true
            x = nodes[i]
            p⃗L = periodic_offset(x_prev, x, Ls, Lhs)
            x = x - p⃗L
            # @assert all(@. abs(x - x_prev) < Lhs)
            push!(xs, x)
            i = node_next[i]
            x_prev = x
        end
        p⃗L = periodic_offset(x_first, x_prev, Ls, Lhs)  # needed to "close" the filament
        Np = length(xs)  # number of points in the filament
        f = Filaments.similar_filament(fbase, (Np,); offset = p⃗L)
        f .= xs
        if !Filaments.check_nodes(Bool, f)
            # Completely remove this filament, as it is too small to be represented by the
            # discretisation method (usually < 3 or < 5 nodes, depending on the method).
            f[end + 1] = f[begin] + p⃗L  # needed for correct length estimation
            local Lfil = Filaments.filament_length(f; quad = nothing)
            filaments_removed_count += 1
            filaments_removed_length += Lfil
        else
            Filaments.update_coefficients!(f)
            push!(fs, f)
            callback(f, lastindex(fs), :appended)
        end
        # Find next node to be reconstructed
        while i_start <= lastindex(reconstructed) && reconstructed[i_start]
            i_start += 1
        end
    end
    (; filaments_removed_count, filaments_removed_length)
end

function _reconnect_pass!(callback::F, ret_base, cache::ReconnectFastCache, fs, vs; to) where {F <: Function}
    (; reconnection_count, reconnection_length_loss, filaments_removed_count, filaments_removed_length,) = ret_base
    rec_before = reconnection_count
    @timeit to "Set points" _update_cache!(cache, fs; to)
    @timeit to "Find reconnections" let
        local ret = _reconnect_from_cache!(cache; to)
        reconnection_count += ret.reconnection_count
        reconnection_length_loss += ret.reconnection_length_loss
    end
    if reconnection_count - rec_before > 0  # if reconnections were performed
        @timeit to "Reconstruct filaments" let
            local ret = _reconstruct_filaments!(callback, fs, cache)
            filaments_removed_count += ret.filaments_removed_count
            filaments_removed_length += ret.filaments_removed_length
        end
    end
    (;
        reconnection_count,
        reconnection_length_loss,
        filaments_removed_count,
        filaments_removed_length,
    )::typeof(ret_base)
end
