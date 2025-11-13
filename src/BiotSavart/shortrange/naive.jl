export NaiveShortRangeBackend

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

struct NaiveShortRangeCache{
        Params <: ParamsShortRange{<:Real, <:NaiveShortRangeBackend},
        Charges <: PointData,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    params :: Params
    data   :: Charges
    to     :: Timer
end

function init_cache_short(
        ::ParamsCommon, params::ParamsShortRange{T, <:NaiveShortRangeBackend},
        data::PointData, to::TimerOutput,
    ) where {T}
    NaiveShortRangeCache(params, data, to)
end

function nearby_charges(c::NaiveShortRangeCache, x⃗::Vec3)
    (; data,) = c
    # Note: it's not worth it to filter out charges that are too far from x⃗, since that job
    # is done again in `biot_savart_contribution`.
    # So we simply return all charges one by one, regardless of x⃗.
    eachindex(data.points, data.charges, data.segments)
end
