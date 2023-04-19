export NaiveShortRangeBackend

"""
    NaiveShortRangeBackend <: ShortRangeBackend

Naive computation of short-range interactions.

This backend can be quite slow and should be used for testing only.
"""
struct NaiveShortRangeBackend <: ShortRangeBackend end

struct NaiveShortRangeCache{
        Params <: ParamsShortRange,
    } <: ShortRangeCache
    params :: Params
end

function init_cache_short(
        common::ParamsCommon{T}, params::ParamsShortRange{<:NaiveShortRangeBackend},
    ) where {T}
    NaiveShortRangeCache(params)
end
