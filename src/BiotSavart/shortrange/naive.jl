export NaiveShortRangeBackend

"""
    NaiveShortRangeBackend <: ShortRangeBackend

Naive computation of short-range interactions.
"""
struct NaiveShortRangeBackend <: ShortRangeBackend end

struct NaiveShortRangeCache{
        Params <: ParamsShortRange{<:NaiveShortRangeBackend},
        Charges <: PointData,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    params :: Params
    data   :: Charges
    to     :: Timer
end

function init_cache_short(
        ::ParamsCommon, params::ParamsShortRange{<:NaiveShortRangeBackend},
        data::PointData, to::TimerOutput,
    )
    NaiveShortRangeCache(params, data, to)
end

function nearby_charges(c::NaiveShortRangeCache, x⃗::Vec3)
    (; data,) = c
    # Note: it's not worth it to filter out charges that are too far from x⃗, since that job
    # is done again in `biot_savart_contribution`.
    # So we simply return all charges one by one, regardless of x⃗.
    Iterators.zip(data.points, data.charges)
end
