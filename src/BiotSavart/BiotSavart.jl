"""
    BiotSavart

Module for estimation of Biot–Savart integrals along vortex filaments using
fast Ewald splitting.
"""
module BiotSavart

export
    ParamsBiotSavart,
    GaussLegendreQuadrature,
    Zero, Infinity, ∞,
    init_cache,
    periods,
    velocity_on_nodes!

using ..BasicTypes:
    Vec3, Derivative, Zero, Infinity, ∞

using ..Quadratures:
    quadrature, GaussLegendreQuadrature, AbstractQuadrature

using ..Filaments:
    Filaments, AbstractFilament, ClosedFilament, CurvatureBinormal,
    knots, nodes, segments, integrate

using TimerOutputs: TimerOutput, @timeit

# Common parameters to short- and long-range computations.
struct ParamsCommon{T, Alpha <: Real, Sigma <: Real, Periods <: NTuple{3, Real}}
    Γ  :: T        # vortex circulation
    a  :: T        # vortex core size
    Δ  :: T        # LIA coefficient given by core vorticity profile
    α  :: Alpha    # Ewald splitting parameter (inverse length scale)
    σ  :: Sigma    # Ewald splitting length scale = 1 / α√2 = std of Gaussian filter
    Ls :: Periods  # size of unit cell (= period in each direction)
    function ParamsCommon{T}(Γ, a, Δ, α, Ls) where {T}
        σ = 1 / (α * sqrt(2))
        new{T, typeof(α), typeof(σ), typeof(Ls)}(Γ, a, Δ, α, σ, Ls)
    end
end

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const VectorOfPositions = AbstractVector{<:Vec3}
const VectorOfVelocities = AbstractVector{<:Vec3}
const AllFilamentVelocities = AbstractVector{<:VectorOfVelocities}

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

- `Γ::Real`: vortex circulation (assumed constant);

- `a::Real`: vortex core size (assumed constant);

- `α::Real`: Ewald splitting parameter (inverse length scale). One can set
  `α = Zero()` to efficiently disable long-range computations.

- `Ls::Union{Real, NTuple{3, Real}}`: size of unit cell (i.e. period in each direction).
  If a single value is passed (e.g. `Ls = 2π`), it is assumed that periods are
  the same in each direction.

  One can set `Ls = ∞` to disable periodicity. This should be done in combination with `α = Zero()`.

- `Ns::Dims{3}`: dimensions of physical grid used for long-range interactions. This parameter
  is not required if `α = Zero()`.

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
            a::Real,
            quadrature_short::AbstractQuadrature = GaussLegendreQuadrature(4),
            quadrature_long::AbstractQuadrature = GaussLegendreQuadrature(2),
            backend_short::ShortRangeBackend = NaiveShortRangeBackend(),
            backend_long::LongRangeBackend = FINUFFTBackend(),
            Δ::Real = 0.25,
            kws...,
        ) where {T}
        # TODO better split into physical (Γ, a, Δ, Ls) and numerical (α, rcut, Ns, ...) parameters?
        # - define ParamsPhysical instead of ParamsCommon
        # - include α in both ParamsShortRange and ParamsLongRange?
        (; Ns, rcut,) = _extra_params(α; kws...)
        common = ParamsCommon{T}(Γ, a, Δ, α, Ls)
        sr = ParamsShortRange(backend_short, quadrature_short, common, rcut)
        lr = ParamsLongRange(backend_long, quadrature_long, common, Ns)
        new{typeof(common), typeof(sr), typeof(lr)}(common, sr, lr)
    end
end

periods(p::ParamsBiotSavart) = p.common.Ls

_extra_params(α::Zero; Ns = (0, 0, 0), rcut = ∞) = (; Ns, rcut,)
_extra_params(α::Real; Ns, rcut = 4 / α) = (; Ns, rcut,)  # Ns is required in this case

ParamsBiotSavart(::Type{T}; Γ::Real, α::Real, Ls, kws...) where {T} =
    ParamsBiotSavart(T, Γ, α, _convert_periods(Ls); kws...)

_convert_periods(Ls::NTuple{3, Real}) = Ls
_convert_periods(L::Real) = (L, L, L)

ParamsBiotSavart(; kws...) = ParamsBiotSavart(Float64; kws...)

"""
    BiotSavartCache

Includes arrays and data required for computation of Biot–Savart integrals.

## Fields

- `shortrange` cache associated to short-range computations;
- `longrange` cache associated to long-range computations. It can be `NullLongRangeCache()`
  in case the Ewald parameter `α` was set to `Zero()`;
- `to` a `TimerOutput` instance for measuring the time spent on different functions.

"""
struct BiotSavartCache{
        ShortRange <: ShortRangeCache,
        LongRange <: LongRangeCache,
        Timer,
    }
    shortrange :: ShortRange
    longrange  :: LongRange
    to         :: Timer
end

"""
    init_cache(p::ParamsBiotSavart; timer = TimerOutput("BiotSavart")) -> BiotSavartCache

Initialise caches for computing Biot–Savart integrals.
"""
function init_cache(p::ParamsBiotSavart; timer = TimerOutput("BiotSavart"))
    shortrange = init_cache_short(p.common, p.shortrange, timer)
    longrange = init_cache_long(p.common, p.longrange, timer)
    BiotSavartCache(shortrange, longrange, timer)
end

"""
    velocity_on_nodes!(
        vs::AbstractVector{<:AbstractVector{<:Vec3}},
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament},
    )

Compute velocity induced by vortex filaments on filament nodes.

Velocities induced by vortex filaments `fs` are written to `vs`.

Usually, `fs` is a vector containing all the vortex filaments in the system.
In that case, `vs` must be a vector of vectors, which will contain the velocities of
all filament nodes. The length of `vs[i]` must be equal to the number of nodes
in the filament `fs[i]`.

For convenience, if the system contains a single vortex filament, one can also
pass a single velocity vector `v` and a single filament `f`.
"""
function velocity_on_nodes! end

function _reset_vectors!(vs)
    for v ∈ vs
        fill!(v, zero(eltype(v)))
    end
    vs
end

function velocity_on_nodes!(
        vs::AllFilamentVelocities,
        cache::BiotSavartCache,
        fs::VectorOfFilaments,
    )
    (; to,) = cache
    eachindex(vs) == eachindex(fs) || throw(DimensionMismatch("wrong dimensions of velocity vector"))
    _reset_vectors!(vs)
    if cache.longrange !== NullLongRangeCache()
        @timeit to "add_long_range_velocity!" add_long_range_velocity!(vs, cache.longrange, fs)
    end
    inds = eachindex(fs)
    @inbounds for (i, f) ∈ pairs(fs)
        for j ∈ first(inds):(i - 1)
            Xs = nodes(fs[j])
            @timeit to "add_short_range_velocity_other!" add_short_range_velocity_other!(vs[j], Xs, cache.shortrange, f)
        end
        @timeit to "add_short_range_velocity_self!" add_short_range_velocity_self!(vs[i], cache.shortrange, f)
        for j ∈ (i + 1):last(inds)
            Xs = nodes(fs[j])
            @timeit to "add_short_range_velocity_other!" add_short_range_velocity_other!(vs[j], Xs, cache.shortrange, f)
        end
    end
    vs
end

# Case of a single filament: interpret inputs as single-element vectors.
function velocity_on_nodes!(
        v::VectorOfVelocities, cache::BiotSavartCache, f::AbstractFilament,
    )
    vs = SVector{1}((v,))
    fs = SVector{1}((f,))
    velocity_on_nodes!(vs, cache, fs)
    v
end

end
