using ..BiotSavart: PointVortexState, Vec2

export AnnihilationCriterion, AbstractAnnihilationCache, NoAnnihilations, AnnihilateFast

"""
    AnnihilationCriterion

Abstract type describing a criterion for annihilation of point vortex dipoles.

Implemented annihilation criteria include:

- [`NoAnnihilations`](@ref): disables annihilations;

- [`AnnihilateFast`](@ref): annihilates vortex dipoles which are closer than a prescribed
  distance.
"""
abstract type AnnihilationCriterion end

abstract type AbstractAnnihilationCache end

function init_cache(crit::AnnihilationCriterion, vortices::PointVortexState, Ls::NTuple{2, Real})
    _init_cache(crit, vortices, Ls)
end

"""
    annihilate!(cache::AbstractAnnihilationCache, vortices::PointVortexState; [folded = false], [to = TimerOutput()]) -> Int

Annihilate point vortex dipoles according to chosen criterion.

This function will remove annihilated vortices from the `vortices` vector.

If all vortex points are within the main periodic cell `[0, L]²`, passing `folded = true`
may slightly accelerate the computations.

Returns the number of removed vortices (`= 2 * n_annihilations`).
"""
function annihilate!(
        cache::AbstractAnnihilationCache,
        vortices::PointVortexState;
        folded = false,
        to::TimerOutput = TimerOutput(),
    )
    _annihilate!(cache, vortices, folded, to)::Int
end

## ==================================================================================================== ##

"""
    NoAnnihilations <: AnnihilationCriterion
    NoAnnihilations()

Used to disable point vortex annihilation.
"""
struct NoAnnihilations <: AnnihilationCriterion end

struct NullAnnihilationCache <: AbstractAnnihilationCache end

_init_cache(::NoAnnihilations, args...) = NullAnnihilationCache()
_annihilate!(::NullAnnihilationCache, args...) = 0

## ==================================================================================================== ##

"""
    AnnihilateFast <: AnnihilationCriterion
    AnnihilateFast(d_crit)

Annihilates point vortex dipoles (with opposite circulations) which are at a distance `d < d_crit`.

Similarly to [`ReconnectFast`](@ref) in 3D, this criterion only supports periodic domains.
"""
struct AnnihilateFast <: AnnihilationCriterion
    dist::Float64
end

distance(c::AnnihilateFast) = c.dist

struct AnnihilateFastCache{
        T <: AbstractFloat,
        Criterion <: AnnihilateFast,
        Finder <: PeriodicCellList,
        IndexVector <: AbstractVector{<:Integer},
    } <: AbstractAnnihilationCache
    crit::Criterion
    cl::Finder
    Ls::NTuple{2, T}
    idx_other::IndexVector  # length Np; index of "other" vortex to annihilate with (0 means no annihilation)
end

criterion(c::AnnihilateFastCache) = c.crit

function Base.show(io::IO, c::AnnihilateFast)
    print(io, "AnnihilateFast($(c.dist))")
end

function _init_cache(crit::AnnihilateFast, vortices, Ls_in)
    (; points) = vortices
    T = eltype(eltype(points))
    Ls = map(T, Ls_in)
    idx_other = Vector{Int}(undef, 0)
    d_annihilate = distance(crit)
    r_cut = T(d_annihilate)
    Lmin = min(Ls...)
    M = 2  # number of subdivisions
    r_cut_max = CellLists.max_cutoff_distance(M, Lmin)
    r_cut > r_cut_max && error(
        lazy"""annihilation distance is too large compared to the domain size: d_annihilate / L_min = $(d_annihilate / Lmin).
        Try a smaller reconnection distance or a larger domain."""
    )
    rs = map(_ -> r_cut, Ls)
    nsubdiv = Val(M)
    cl = PeriodicCellList(rs, Ls, nsubdiv)  # a single element is a point vortex index (in 1:Np)
    AnnihilateFastCache(crit, cl, Ls, idx_other)
end

function _update_cache!(cache::AnnihilateFastCache, points, folded::Bool)
    (; cl, idx_other) = cache
    Base.require_one_based_indexing(points)
    resize!(idx_other, length(points))
    Threads.@threads for i in eachindex(idx_other)
        @inbounds idx_other[i] = 0
    end
    if folded
        CellLists.set_elements!(cl, points; folded = Val(true))  # CellLists can assume that points are already in [0, L]²
    else
        CellLists.set_elements!(cl, points; folded = Val(false))
    end
    cache
end

function _find_nearby_dipoles!(cache, vortices, folded::Val)
    (; cl, crit, idx_other) = cache
    (; points, signs) = vortices
    r_cut = distance(crit)
    r²_cut = r_cut * r_cut
    CellLists.foreach_pair(cl, points; folded = Val(true)) do x⃗, i, j
        @inline
        @inbounds y⃗ = points[j]
        is_dipole = @inbounds iszero(signs[i] + signs[j])  # opposite signs
        r² = sum(abs2, x⃗ - y⃗)
        if is_dipole && r² < r²_cut
            # Note: we don't have the right of modifying idx_other[j] in parallel without atomics.
            # We will let a different call do it, e.g. when the same pair is visited in
            # reverse order (i, j) <-> (j, i).
            n = @inbounds idx_other[i]
            if n == 0
                # No (i, j) pair had been found previously
                @inbounds idx_other[i] = j
            else
                # Check if this pair is "better" (closer) than that found previously.
                z⃗ = @inbounds points[n]
                r²_prev = sum(abs2, x⃗ - z⃗)
                if r²_prev > r²
                    @inbounds idx_other[i] = j
                end
            end
        end
    end
    nothing
end

function _annihilate_dipoles!(cache, vortices)
    (; idx_other) = cache
    (; points, signs) = vortices

    # First process the list of pairs to skip "asymmetric" annihilations (when vortex j is
    # the closest to i, but i is not the closest to j).
    # In that case we set idx_other[i] = 0, which means it won't be removed.
    n_annihilated = 0
    @inbounds for i in eachindex(idx_other)
        j = idx_other[i]
        # Only annihilate if i == idx_other[idx_other[i]], i.e. if both vortices "agree"
        # that this is the best pair. Otherwise we simply skip the annihilation.
        if i != idx_other[j]
            idx_other[i] = 0  # don't remove vortex i
        elseif i < j
            # Annihilate vortices i and j
            n_annihilated += 2
        end
    end

    # Now remove vortices.
    n_annihilated == 0 && return n_annihilated

    isave = firstindex(points) - 1
    @inbounds for n in eachindex(points, signs, idx_other)
        keep = idx_other[n] == 0
        if keep
            isave += 1
            points[isave] = points[n]
            signs[isave] = signs[n]
        end
    end
    @assert isave == length(points) - n_annihilated
    resize!(points, isave)
    resize!(signs, isave)

    n_annihilated
end

function _annihilate!(cache::AnnihilateFastCache, vortices, folded::Bool, to)
    @timeit to "Set points" _update_cache!(cache, vortices.points, folded)
    @timeit to "Find nearby dipoles" if folded
        _find_nearby_dipoles!(cache, vortices, Val(true))
    else
        _find_nearby_dipoles!(cache, vortices, Val(false))
    end
    @timeit to "Annihilate" n_annihilated = _annihilate_dipoles!(cache, vortices)
    n_annihilated
end
