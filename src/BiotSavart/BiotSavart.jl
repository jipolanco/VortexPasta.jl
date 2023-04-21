"""
    BiotSavart

Module for estimation of Biot–Savart integrals along vortex filaments using
fast Ewald splitting.
"""
module BiotSavart

export
    ParamsBiotSavart,
    GaussLegendreQuadrature

using ..BasicTypes:
    Vec3, Derivative

using ..Quadratures:
    quadrature, GaussLegendreQuadrature, AbstractQuadrature

using ..Filaments:
    Filaments, AbstractFilament, ClosedFilament, knots, segments

# Common parameters to short- and long-range computations.
struct ParamsCommon{T}
    Γ  :: T             # vortex circulation
    a  :: T             # vortex core size
    Δ  :: T             # LIA coefficient given by core vorticity profile
    α  :: T             # Ewald splitting parameter (inverse length scale)
    σ  :: T             # Ewald splitting length scale = 1 / α√2 = std of Gaussian filter
    Ls :: NTuple{3, T}  # size of unit cell (= period in each direction)
    function ParamsCommon{T}(Γ, a, Δ, α, Ls) where {T}
        σ = 1 / (α * sqrt(2))
        new{T}(Γ, a, Δ, α, σ, Ls)
    end
end

include("shortrange/shortrange.jl")
include("longrange/longrange.jl")

"""
    ParamsBiotSavart{T <: AbstractFloat}

Contains parameters for calculation of Biot–Savart integrals using fast Ewald splitting.

The type parameter `T` corresponds to the precision used in computations
(typically `Float64` or `Float32`).

# Construction

    ParamsBiotSavart([T = Float64]; Γ, α, Ls, Ns, optional_kws...)

where the optional parameter `T` sets the numerical precision.

Mandatory and optional keyword arguments are detailed in the following.

## Mandatory keyword arguments

- `Γ::Real` vortex circulation (assumed constant);

- `a::Real` vortex core size (assumed constant);

- `α::Real` Ewald splitting parameter (inverse length scale);

- `Ls::NTuple{3, Real}` size of unit cell (i.e. period in each direction);

- `Ns::Dims{3}` dimensions of physical grid used for long-range interactions.

## Optional keyword arguments (and their defaults)

### Short-range interactions

- `backend_short::ShortRangeBackend = NaiveShortRangeBackend()` backend used to compute
  short-range interactions;

- `quadrature_short::AbstractQuadrature = GaussLegendreQuadrature(4)`
  quadrature rule for short-range interactions;

- `rcut = 4√2 / α` cutoff distance for computation of short-range interactions.
  For performance reasons, the cutoff distance must be less than half the cell
  unit size in each direction, i.e. `rcut < minimum(Ls) / 2`.

### Long-range interactions

- `backend_long::LongRangeBackend = FINUFFTBackend()` backend used to compute
  long-range interactions;

- `quadrature_long::AbstractQuadrature = GaussLegendreQuadrature(2)`
  quadrature rule for long-range interactions.

### Local self-induced velocity

- `Δ = 0.25` coefficient appearing in the local self-induced velocity (LIA
  term), which depends on the vorticity profile at the vortex core.

  Some common values of `Δ` are:

  * `Δ = 0.25` for a constant vorticity profile (default);

  * `Δ = 0.5` for a hollow vortex;

  * `Δ ≈ 0.905 ≈ 0.558 + ln(2) / 2` for a Gaussian vorticity profile (taking
    `a` as the Gaussian standard deviation `σ`);

  * `Δ ≈ 0.615` for a Gross–Pitaevskii vortex with healing length `a`.

  See Saffman (1992), sections 10.2--10.3 for the first three.

"""
struct ParamsBiotSavart{
        Common <: ParamsCommon,
        ShortRange <: ParamsShortRange,
        LongRange <: ParamsLongRange,
    }
    common     :: Common
    shortrange :: ShortRange
    longrange  :: LongRange

    function ParamsBiotSavart(
            ::Type{T}, Γ::Real, α::Real, Ls::NTuple{3, Real};
            a::Real, Ns::Dims{3}, 
            quadrature_short::AbstractQuadrature = GaussLegendreQuadrature(4),
            quadrature_long::AbstractQuadrature = GaussLegendreQuadrature(2),
            backend_short::ShortRangeBackend = NaiveShortRangeBackend(),
            backend_long::LongRangeBackend = FINUFFTBackend(),
            Δ::Real = 0.25,
            rcut = 4 / α,
        ) where {T}
        # TODO better split into physical (Γ, a, Δ, Ls) and numerical (α, rcut, Ns, ...) parameters?
        # - define ParamsPhysical instead of ParamsCommon
        # - include α in both ParamsShortRange and ParamsLongRange?
        common = ParamsCommon{T}(Γ, a, Δ, α, Ls)
        sr = ParamsShortRange(backend_short, quadrature_short, common, rcut)
        lr = ParamsLongRange(backend_long, quadrature_long, Ns)
        new{typeof(common), typeof(sr), typeof(lr)}(common, sr, lr)
    end
end

ParamsBiotSavart(::Type{T}; Γ::Real, α::Real, Ls::NTuple, kws...) where {T} =
    ParamsBiotSavart(T, Γ, α, Ls; kws...)

ParamsBiotSavart(; kws...) = ParamsBiotSavart(Float64; kws...)

"""
    BiotSavartCache

Includes arrays and data required for computation of Biot–Savart integrals.

## Fields

- `longrange` cache associated to long-range computations
"""
struct BiotSavartCache{
        ShortRange <: ShortRangeCache,
        LongRange <: LongRangeCache,
    }
    shortrange :: ShortRange
    longrange  :: LongRange
end

"""
    init_cache(p::ParamsBiotSavart) -> BiotSavartCache

Initialise caches for computing Biot–Savart integrals.
"""
function init_cache(p::ParamsBiotSavart)
    shortrange = init_cache_short(p.common, p.shortrange)
    longrange = init_cache_long(p.common, p.longrange)
    BiotSavartCache(shortrange, longrange)
end

end
