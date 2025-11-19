struct ShortRangeCacheCommon{
        Params <: ParamsShortRange,
        PointCharges <: PointData,
        Timer <: TimerOutput,
    }
    params    :: Params
    pointdata :: PointCharges
    to        :: Timer
end
