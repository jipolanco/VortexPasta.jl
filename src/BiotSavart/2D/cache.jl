"""
    PointVortexCache

Includes arrays and data required for computation of point vortex Biot–Savart interactions.

## Fields

- `pointdata`: contains vector values at points in space. This is used for both short-range
  and long-range computations;

- `shortrange`: cache associated to short-range computations;

- `longrange`: cache associated to long-range computations;

- `to`: a `TimerOutput` instance for measuring the time spent on different functions.

"""
struct PointVortexCache{
        Params <: ParamsPointVortex,
        Data <: PointVortexData,
        ShortRange <: ShortRangeCache,
        LongRange <: LongRangeCache,
        Timer,
    }
    params     :: Params
    pointdata  :: Data
    shortrange :: ShortRange
    longrange  :: LongRange
    to         :: Timer
end

Base.summary(io::IO, c::PointVortexCache) = print(io, "PointVortexCache")

function init_cache(
        p::ParamsPointVortex;  # argument 2 is for backwards compatibility (no longer needed)
        timer = TimerOutput("PointVortex"),
    )
    T = eltype(p)
    pointdata = PointVortexData(T)
    shortrange = init_cache_short(p, pointdata)
    longrange = init_cache_long(p, pointdata)
    PointVortexCache(p, pointdata, shortrange, longrange, timer)
end

# This is called from LongRangeCacheCommon
_longrange_outputs(data::PointVortexData{T}) where {T} = (;
    velocity = similar(data.points),           # 2D vectors
    streamfunction = similar(data.points, T),  # scalars
)

# This is called from ShortRangeCacheCommon
_shortrange_outputs(data::PointVortexData{T}) where {T} = (;
    velocity = similar(data.points),           # 2D vectors
    streamfunction = similar(data.points, T),  # scalars
)
