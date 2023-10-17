export NaiveShortRangeBackend

using ..FindNearbySegments: NaiveSegmentFinder

"""
    NaiveShortRangeBackend <: ShortRangeBackend

Naive computation of short-range interactions.
"""
struct NaiveShortRangeBackend <: ShortRangeBackend end

struct NaiveShortRangeCache{
        Finder <: NaiveSegmentFinder,
        Params <: ParamsShortRange{<:NaiveShortRangeBackend},
        Timer <: TimerOutput,
    } <: ShortRangeCache
    finder :: Finder
    params :: Params
    to     :: Timer
end

function init_cache_short(
        ::ParamsCommon, params::ParamsShortRange{<:NaiveShortRangeBackend},
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    )
    (; common, rcut,) = params
    (; Ls,) = common
    finder = NaiveSegmentFinder(fs, rcut, Ls)
    NaiveShortRangeCache(finder, params, to)
end
