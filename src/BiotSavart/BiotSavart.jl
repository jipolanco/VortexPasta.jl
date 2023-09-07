"""
    BiotSavart

Module for estimation of Biot–Savart integrals along vortex filaments using
fast Ewald splitting.
"""
module BiotSavart

export
    ParamsBiotSavart,
    GaussLegendre,
    Zero, Infinity, ∞,
    Velocity, Streamfunction,
    init_cache,
    periods,
    velocity_on_nodes!,
    compute_on_nodes!

using ..BasicTypes:
    Vec3, Derivative, Zero, Infinity, ∞

using ..Quadratures:
    quadrature, GaussLegendre, AbstractQuadrature

using ..Filaments:
    Filaments, AbstractFilament, ClosedFilament, CurvatureBinormal,
    knots, nodes, segments, integrate

using TimerOutputs: TimerOutput, @timeit

abstract type OutputField end
struct Streamfunction <: OutputField end
struct Velocity <: OutputField end

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

Base.eltype(::Type{<:ParamsCommon{T}}) where {T} = T
Base.eltype(p::ParamsCommon) = eltype(typeof(p))

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const VectorOfVec = AbstractVector{<:Vec3}
const VectorOfPositions = VectorOfVec
const VectorOfVelocities = VectorOfVec
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

  One can set `Ls = Infinity()` to disable periodicity. This should be done in combination with `α = Zero()`.

- `Ns::Dims{3}`: dimensions of physical grid used for long-range interactions. This parameter
  is not required if `α = Zero()`.

## Optional keyword arguments (and their defaults)

### Short-range interactions

- `backend_short::ShortRangeBackend`: backend used to compute
  short-range interactions. The default is `CellListsBackend(2)`, unless periodicity is
  disabled, in which case `NaiveShortRangeBackend()` is used.
  See [`ShortRangeBackend`](@ref) for a list of possible backends;

- `quadrature_short::AbstractQuadrature = GaussLegendre(4)`:
  quadrature rule for short-range interactions;

- `rcut = 4√2 / α`: cutoff distance for computation of short-range interactions.
  For performance and practical reasons, the cutoff distance must be less than half the cell
  unit size in each direction, i.e. `rcut < minimum(Ls) / 2`.

### Long-range interactions

- `backend_long::LongRangeBackend = FINUFFTBackend()`: backend used to compute
  long-range interactions. See [`LongRangeBackend`](@ref) for a list of possible backends;

- `quadrature_long::AbstractQuadrature = GaussLegendre(2)`:
  quadrature rule for long-range interactions.

### Local self-induced velocity

- `Δ = 0.25`: coefficient appearing in the local self-induced velocity (LIA
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
        T,
        Common <: ParamsCommon{T},
        ShortRange <: ParamsShortRange,
        LongRange <: ParamsLongRange,
    }
    common     :: Common
    shortrange :: ShortRange
    longrange  :: LongRange

    function ParamsBiotSavart(
            ::Type{T}, Γ::Real, α::Real, Ls::NTuple{3, Real};
            a::Real,
            quadrature_short::AbstractQuadrature = GaussLegendre(4),
            quadrature_long::AbstractQuadrature = GaussLegendre(2),
            backend_short::ShortRangeBackend = default_short_range_backend(Ls),
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
        new{T, typeof(common), typeof(sr), typeof(lr)}(common, sr, lr)
    end
end

# Returns `true` if `Ls` contains `Infinity` (one or more times), `false` otherwise.
is_open_domain(Ls::Tuple) = is_open_domain(Ls...)
is_open_domain(::Infinity, etc...) = true
is_open_domain(::Real, etc...) = is_open_domain(etc...)
is_open_domain() = false

function default_short_range_backend(Ls::Tuple)
    if is_open_domain(Ls)
        NaiveShortRangeBackend()
    else
        CellListsBackend(2)  # use 2 subdivisions by default (generally faster)
    end
end

# Returns the float type used (e.g. Float64)
Base.eltype(::Type{<:ParamsBiotSavart{T}}) where {T} = T
Base.eltype(p::ParamsBiotSavart) = eltype(typeof(p))

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

- `shortrange`: cache associated to short-range computations;
- `longrange`: cache associated to long-range computations. It can be `NullLongRangeCache()`
  in case the Ewald parameter `α` was set to `Zero()`;
- `to`: a `TimerOutput` instance for measuring the time spent on different functions.

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
    shortrange = init_cache_short(p.common, p.shortrange, fs, timer)
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

This is the same as calling [`compute_on_nodes!`](@ref) when only the velocity is needed.

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
        vs::AbstractVector{<:VectorOfVec},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
        LIA = Val(true),
    )
    fields = (; velocity = vs,)
    compute_on_nodes!(fields, cache, fs; LIA)
end

# Extracts and resets output fields (velocity and/or streamfunction).
function _setup_fields!(fields::NamedTuple, fs)
    N = length(fields)
    @assert N > 0
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    vecs = filter(!isnothing, (vs, ψs))
    nfields = length(vecs)
    nfields == N || throw(ArgumentError(lazy"some of these fields were not recognised: $Names"))
    for us ∈ vecs
        eachindex(us) == eachindex(fs) || throw(DimensionMismatch("wrong dimensions of vector"))
        _reset_vectors!(us)
    end
    (; vs, ψs,)
end

# Returns a tuple (Streamfunction() => ψs, Velocity() => vs) if both fields are available.
# Otherwise, it returns a subset of this which only includes the available fields.
function _fields_to_pairs(fields::NamedTuple{Names}) where {Names}
    N = length(fields)
    @assert N > 0
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    ps = (
        _make_field_pair(Velocity(), vs)...,
        _make_field_pair(Streamfunction(), ψs)...,
    )
    length(ps) == N || throw(ArgumentError(lazy"some of these fields were not recognised: $Names"))
    ps
end

_make_field_pair(key, ::Nothing) = ()
_make_field_pair(key, vs) = (key => vs,)

"""
    compute_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament};
        LIA = Val(true),
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}

Compute velocity and/or streamfunction on filament nodes.

The first argument contains one or more output fields to compute. It is usually of length 1
or 2, and can contain fields named `velocity` and `streamfunction`.

For example, to compute both velocity and streamfunction on the nodes of filaments `fs`:

```julia
# Initialise fields to compute (vectors of vectors)
vs = map(similar ∘ nodes, fs)  # one velocity vector per filament node
ψs = map(similar, vs)

# The first argument to `compute_on_nodes!` must have the following form.
# One can also choose to pass just one of the two fields.
fields = (;
    velocity = vs,
    streamfunction = ψs,
)

cache = BiotSavart.init_cache(...)
compute_on_nodes!(fields, cache, fs)
```

## Disabling LIA / computing *only* LIA

One may disable computation of the locally-induced velocity and streamfunction (LIA term)
by passing `LIA = Val(false)`. Conversely, one can pass `LIA = Val(:only)` to compute *only*
the LIA term. This can be useful for splitting the induced filament
velocities/streamfunctions onto local and non-local parts.
"""
function compute_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
        LIA = Val(true),
        longrange = Val(true),
        shortrange = Val(true),
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    if LIA === Val(:only)
        return _compute_LIA_on_nodes!(fields, cache, fs)
    end

    (; to,) = cache
    (; vs, ψs,) = _setup_fields!(fields, fs)

    if cache.longrange !== NullLongRangeCache() && longrange === Val(true)
        @timeit to "Long-range component" begin
            @timeit to "Vorticity to Fourier" compute_vorticity_fourier!(cache.longrange, fs)
            set_interpolation_points!(cache.longrange, fs)
            if ψs !== nothing
                @timeit to "Streamfunction" begin
                    to_smoothed_streamfunction!(cache.longrange)
                    interpolate_to_physical!(cache.longrange)
                    add_long_range_output!(ψs, cache.longrange)
                end
            end
            if vs !== nothing
                # Velocity must be computed after streamfunction if both are enabled.
                @timeit to "Velocity" begin
                    to_smoothed_velocity!(cache.longrange)
                    interpolate_to_physical!(cache.longrange)
                    add_long_range_output!(vs, cache.longrange)
                end
            end
        end
    end

    if shortrange === Val(true)
        @timeit to "Short-range component" begin
            set_filaments!(cache.shortrange, fs)
            for i ∈ eachindex(fs)
                fields_i = map(us -> us[i], fields)  # velocity/streamfunction of i-th filament
                add_short_range_fields!(fields_i, cache.shortrange, fs[i]; LIA)
            end
        end
    end

    fields
end

function _compute_LIA_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    (; to,) = cache
    (; params,) = cache.shortrange
    # Note: we must use the same quadrature as used when computing the globally induced terms
    (; quad,) = params
    (; Γ, a, Δ,) = params.common
    ps = _fields_to_pairs(fields)  # e.g. (Velocity() => vs, Streamfunction() => ψs)
    prefactor = Γ / (4π)
    @timeit to "LIA term (only)" begin
        @inbounds for (n, f) ∈ pairs(fs), i ∈ eachindex(f)
            for (quantity, values) ∈ ps
                # Here `quantity` is either Velocity() or Streamfunction()
                values[n][i] = BiotSavart.local_self_induced(
                    quantity, f, i, prefactor;
                    a, Δ, quad,
                )
            end
        end
    end
    fields
end

end
