export NaiveShortRangeBackend
using StaticArrays: MVector

@doc raw"""
    NaiveShortRangeBackend <: ShortRangeBackend

Naive computation of short-range interactions.

## Maximum cut-off distance

In periodic domains, this backend requires a cut-off distance ``r_{\text{cut}}`` not larger
than half the domain period ``L`` in each direction:

```math
r_{\text{cut}} ≤ \frac{L}{2}
```
"""
struct NaiveShortRangeBackend <: ShortRangeBackend end

struct NaiveShortRangeCache{Common <: ShortRangeCacheCommon} <: ShortRangeCache
    common :: Common
end

function init_cache_short(
        ::ParamsCommon, params::ParamsShortRange{T, <:NaiveShortRangeBackend}, pointdata::PointData,
    ) where {T}
    common = ShortRangeCacheCommon(params, pointdata)
    NaiveShortRangeCache(common)
end

@inline function nearby_charges(c::NaiveShortRangeCache, x⃗::Vec3)
    (; pointdata,) = c
    # Note: it's not worth it to filter out charges that are too far from x⃗, since that job
    # is done again in `biot_savart_contribution`.
    # So we simply return all charges one by one, regardless of x⃗.
    eachindex(pointdata.points, pointdata.charges, pointdata.segments)
end

@inline function foreach_charge(
        f::F, c::NaiveShortRangeCache, x⃗::Vec3;
        batch_size::Union{Nothing, Val} = nothing,  # see CellLists.foreach_source for details on this argument
        folded = Val(false),  # ignored (used in CellLists)
    ) where {F <: Function}
    _foreach_charge(f, batch_size, c, x⃗)
    nothing
end

@inline function _foreach_charge(
        f::F, batch_size::Nothing, c::NaiveShortRangeCache, x⃗::Vec3
    ) where {F}
    for j in nearby_charges(c, x⃗)
        @inline f(j)
    end
    nothing
end

@inline function _foreach_charge(
        f::F, ::Val{batch_size}, c::NaiveShortRangeCache, x⃗::Vec3
    ) where {F, batch_size}
    inds = MVector{batch_size, Int}(undef)
    m = 0
    for j in nearby_charges(c, x⃗)
        @inbounds inds[m += 1] = j
        if m == batch_size
            @inline f(Tuple(inds), batch_size)
            m = 0
        end
    end
    if m > 0
        for l in (m + 1):batch_size
            # copy latest index, just to make sure that all returned indices are valid
            @inbounds inds[l] = inds[m]
        end
        @inline f(Tuple(inds), m)
    end
    nothing
end

@inline function foreach_pair(
        f::F, c::NaiveShortRangeCache;
        batch_size::Union{Nothing, Val} = nothing,  # see CellLists.foreach_source for details on this argument
        folded = Val(false),  # ignored (used in CellLists)
    ) where {F <: Function}
    (; nodes,) = c.pointdata
    Threads.@threads for i in eachindex(nodes)
        x⃗ = @inbounds nodes[i]
        g = CellLists.Fix12(f, x⃗, i)
        _foreach_charge(g, batch_size, c, x⃗)
    end
    nothing
end
