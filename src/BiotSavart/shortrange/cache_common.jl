struct ShortRangeCacheCommon{
        Params <: ParamsShortRange,
        Splitting <: AbstractEwaldSplitting,
        PointCharges <: PointData,
        OutputVectors <: NamedTuple,
        Timer <: TimerOutput,
    }
    params    :: Params
    splitting :: Splitting      # copy of params.common.splitting, but with arrays on the device (e.g. for Chebyshev series evaluations, if KaiserBesselSplitting)
    pointdata :: PointCharges   # points, charges and evaluation nodes (filament discretisation points)
    outputs   :: OutputVectors  # output velocity and streamfunction fields (as linear vectors)
    to        :: Timer
end

function ShortRangeCacheCommon(params::ParamsShortRange, pointdata_in::PointData)
    (; backend,) = params
    (; splitting,) = params.common
    ka_backend = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    # Make sure we've activated the device (e.g. GPU id) where short-range computations will
    # be performed. We need arrays to be allocated in that device.
    expected_device = KA.device(backend)  # 1, 2, ...
    @assert KA.device(ka_backend) == expected_device
    pointdata = adapt(ka_backend, pointdata_in)      # create PointData replica on the device if needed
    pointdata = if ka_backend isa CPU
        copy(pointdata_in)  # make sure pointdata and pointdata_in are not aliased!
    else
        adapt(ka_backend, pointdata_in)  # create PointData replica on the device
    end
    splitting_d = adapt(ka_backend, splitting)  # splitting
    outputs = (;
        velocity = similar(pointdata.charges),
        streamfunction = similar(pointdata.charges),
    )
    to = TimerOutput()
    ShortRangeCacheCommon(params, splitting_d, pointdata, outputs, to)
end
