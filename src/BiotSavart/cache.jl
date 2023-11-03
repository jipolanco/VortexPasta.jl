"""
    BiotSavartCache

Includes arrays and data required for computation of Biot–Savart integrals.

## Fields

- `pointdata`: contains vector values at points in space. This is used for both short-range
  and long-range computations;

- `shortrange`: cache associated to short-range computations;

- `longrange`: cache associated to long-range computations. It can be `NullLongRangeCache()`
  in case the Ewald parameter `α` was set to `Zero()`;

- `to`: a `TimerOutput` instance for measuring the time spent on different functions.

"""
struct BiotSavartCache{
        Params <: ParamsBiotSavart,
        Data <: PointData,
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

Base.summary(io::IO, c::BiotSavartCache) = print(io, "BiotSavartCache")

"""
    init_cache(
        p::ParamsBiotSavart, fs::AbstractVector{<:AbstractFilament};
        timer = TimerOutput("BiotSavart"),
    ) -> BiotSavartCache

Initialise caches for computing Biot–Savart integrals.
"""
function init_cache(
        p::ParamsBiotSavart, fs::AbstractVector{<:AbstractFilament};
        timer = TimerOutput("BiotSavart"),
    )
    T = eltype(p)
    pointdata = PointData(T)
    shortrange = init_cache_short(p.common, p.shortrange, pointdata, timer)
    longrange = init_cache_long(p.common, p.longrange, pointdata, timer)
    BiotSavartCache(p, pointdata, shortrange, longrange, timer)
end
