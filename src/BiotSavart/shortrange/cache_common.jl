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

function ShortRangeCacheCommon(params::ParamsShortRange, pointdata_in::PointData)
    (; backend,) = params
    ka_backend = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    # Make sure we've activated the device (e.g. GPU id) where short-range computations will
    # be performed. We need arrays to be allocated in that device.
    expected_device = KA.device(backend)  # 1, 2, ...
    @assert KA.device(ka_backend) == expected_device
    pointdata = adapt(ka_backend, pointdata_in)      # create PointData replica on the device if needed
    if pointdata === pointdata_in       # basically if ka_backend isa CPU
        pointdata = copy(pointdata_in)  # make sure pointdata and pointdata_in are not aliased!
    end
    outputs = (;
        velocity = similar(pointdata.charges),
        streamfunction = similar(pointdata.charges),
    )
    to = TimerOutput()
    ShortRangeCacheCommon(params, pointdata, outputs, to)
end
