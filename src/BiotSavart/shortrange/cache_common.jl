struct ShortRangeCacheCommon{
        Params <: ParamsShortRange,
        PointCharges <: PointData,
        OutputVectors <: NamedTuple,
        Timer <: TimerOutput,
    }
    params    :: Params
    pointdata :: PointCharges   # points, charges and evaluation nodes (filament discretisation points)
    outputs   :: OutputVectors  # output velocity and streamfunction fields (as linear vectors)
    to        :: Timer
end

function ShortRangeCacheCommon(params::ParamsShortRange, pointdata_in::PointData, to)
    (; backend,) = params
    ka_backend = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    pointdata = adapt(ka_backend, pointdata_in)      # create PointData replica on the device if needed
    if pointdata === pointdata_in       # basically if ka_backend isa CPU
        pointdata = copy(pointdata_in)  # make sure pointdata and pointdata_in are not aliased!
    end
    outputs = (;
        velocity = similar(pointdata.charges),
        streamfunction = similar(pointdata.charges),
    )
    ShortRangeCacheCommon(params, pointdata, outputs, to)
end
