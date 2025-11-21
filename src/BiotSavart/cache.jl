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
    init_cache(p::ParamsBiotSavart; timer = TimerOutput("BiotSavart")) -> BiotSavartCache

Initialise caches for computing Biot–Savart integrals.
"""
function init_cache(
        p::ParamsBiotSavart, fs = nothing;  # argument 2 is for backwards compatibility (no longer needed)
        timer = TimerOutput("BiotSavart"),
    )
    T = eltype(p)
    pointdata = PointData(T)
    shortrange = init_cache_short(p.common, p.shortrange, pointdata)
    longrange = init_cache_long(p, pointdata)
    BiotSavartCache(p, pointdata, shortrange, longrange, timer)
end

"""
    BiotSavart.get_longrange_field_fourier(cache) -> NamedTuple

Obtain long-range field used to compute long-range interactions.

The input `cache` can be a [`BiotSavartCache`](@ref) or a [`LongRangeCache`](@ref).

This function returns a `NamedTuple` with the fields:

- `field`: a `Tuple` `(ux, uy, uz)` describing a vector field in Fourier space.
  Each component is a 3D matrix of complex values.

- `wavenumbers`: a `Tuple` `(kx, ky, kz)` with the wavenumbers associated to the grid in
  Fourier space.

- `state`: allows to know what the returned field actually represents (vorticity, velocity, ...).
  See [`LongRangeCacheState`](@ref) for details.
"""
function get_longrange_field_fourier end

get_longrange_field_fourier(cache::BiotSavartCache) = get_longrange_field_fourier(cache.longrange)

function get_longrange_field_fourier(longrange::LongRangeCache)
    (; uhat, wavenumbers, state,) = longrange.common
    uhat_tup = StructArrays.components(uhat)::NTuple
    (; field = uhat_tup, wavenumbers = wavenumbers, state = copy(state),)
end
