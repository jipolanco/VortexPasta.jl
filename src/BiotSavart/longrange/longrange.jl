using StaticArrays: SVector
using StructArrays: StructArrays, StructVector, StructArray

## ================================================================================ ##

"""
    LongRangeBackend

Abstract type denoting the backend to use for computing long-range interactions.

# Implemented backends

- [`FINUFFTBackend`](@ref): estimates long-range interactions via the non-uniform fast
  Fourier transform (NUFFT) using the FINUFFT library;

- [`ExactSumBackend`](@ref): computes long-range interactions using exact Fourier sums. This
  is really inefficient and should only be used for testing.

# Extended help

## Implementation details

The following functions must be implemented by a `BACKEND <: LongRangeBackend`:

- `init_cache_long_ewald(c::ParamsCommon, p::ParamsLongRange{<:BACKEND}, to::TimerOutput) -> LongRangeCache`.

- [`expected_period`](@ref) (optional),

- [`folding_limits`](@ref) (optional).

"""
abstract type LongRangeBackend end

"""
    expected_period(::LongRangeBackend) -> Union{Nothing, Real}

Domain period expected by the backend.

This is used for rescaling input coordinates to the requirements of the backend.
For instance, FINUFFT assumes a period ``2π``, and therefore coordinates are
rescaled if the input data has a period different from ``2π``.
"""
expected_period(::LongRangeBackend) = nothing

"""
    folding_limits(::LongRangeBackend) -> Union{Nothing, NTuple{2, Real}}

Domain limits required by the backend.

This is used for folding input coordinates so that they are within the limits
expected by the backend.
For instance, FINUFFT requires coordinates to be in the ``[-3π, 3π]`` interval.

Note that, if a backend defines `folding_limits`, then it must also define
[`expected_period`](@ref).
"""
folding_limits(::LongRangeBackend) = nothing

## ================================================================================ ##

struct ParamsLongRange{
        Backend <: LongRangeBackend,
        Quadrature <: AbstractQuadrature,
        Common <: ParamsCommon,
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    Ns      :: Dims{3}     # grid dimensions for FFTs
end

backend(p::ParamsLongRange) = p.backend
quadrature_rule(p::ParamsLongRange) = p.quad

## ================================================================================ ##

mutable struct LongRangeCacheState
    quantity :: Symbol  # quantity currently held by the cache (:undef, :vorticity, :velocity, :streamfunction)
    smoothed :: Bool    # true if Ewald's Gaussian filter has already been applied
end

LongRangeCacheState() = LongRangeCacheState(:undef, false)

struct LongRangeCacheCommon{
        T <: AbstractFloat,
        Params <: ParamsLongRange,
        WaveNumbers <: NTuple{3, AbstractVector},
        Charges <: StructVector{Vec3{Complex{T}}},
        FourierVectorField <: StructArray{Vec3{Complex{T}}, 3},
        Timer <: TimerOutput,
    }
    params      :: Params
    wavenumbers :: WaveNumbers
    charges     :: Charges  # values at non-uniform locations (3 × [Np]) -- only used to store interpolations
    uhat        :: FourierVectorField  # uniform Fourier-space data (3 × [Nx, Ny, Nz])
    ewald_op    :: Array{T, 3}  # Ewald operator in Fourier space ([Nx, Ny, Nz])
    state       :: LongRangeCacheState
    to          :: Timer
end

function LongRangeCacheCommon(
        pcommon::ParamsCommon,
        params::ParamsLongRange,
        wavenumbers::NTuple{3, AbstractVector},
        timer::TimerOutput,
    )
    T = eltype(pcommon)
    (; Γ, Ls, α,) = pcommon
    @assert α !== Zero()
    Nks = map(length, wavenumbers)
    charges = StructVector{Vec3{Complex{T}}}(undef, 0)
    uhat = StructArray{Vec3{Complex{T}}}(undef, Nks)
    ewald_op = init_ewald_fourier_operator(T, wavenumbers, Γ, α, Ls)
    state = LongRangeCacheState()
    LongRangeCacheCommon(params, wavenumbers, charges, uhat, ewald_op, state, timer)
end

"""
    LongRangeCache

Abstract type describing the storage of data required to compute long-range interactions.

The [`init_cache_long`](@ref) function returns a concrete instance of a `LongRangeCache`
(or `NullLongRangeCache()`, if long-range computations were disabled by setting `α = Zero()`).

# Extended help

## Implementation details

### Fields

All caches must include a `common <: LongRangeCacheCommon` field which contains common
definitions for all backends.

Optional fields:

- `points :: StructVector{Vec3{T}}` contains non-uniform points. If the backend defines
  [`expected_period`](@ref), then the `points` field *must* be defined.

### Functions

The following functions must be implemented by a cache:

- [`reset_fields!`](@ref) (optional),

- [`set_num_points!`](@ref),

- [`add_point!`](@ref),

- [`add_pointcharge!`](@ref),

- [`transform_to_fourier!`](@ref),

- [`interpolate_to_physical!`](@ref).
"""
abstract type LongRangeCache end

# This is for convenience, to treat the fields of c.common as if they were fields of c.
# Allows to do things like c.wavenumbers as a shortcut for c.common.wavenumbers.
Base.propertynames(c::LongRangeCache) =
    (fieldnames(typeof(c))..., propertynames(c.common)...)

Base.getproperty(c::LongRangeCache, name::Symbol) = _getproperty(c, Val(name))

function _getproperty(c::LongRangeCache, ::Val{name}) where {name}
    if hasfield(typeof(c), name)
        getfield(c, name)
    else
        getproperty(c.common, name)
    end
end

"""
    NullLongRangeCache <: LongRangeCache

Dummy cache type returned by [`init_cache_long`](@ref) when long-range
computations are disabled.

This is the case when the Ewald splitting parameter ``α`` is set to `Zero()`.
"""
struct NullLongRangeCache <: LongRangeCache end

backend(c::LongRangeCache) = backend(c.params)

"""
    init_cache_long(pc::ParamsCommon, p::ParamsLongRange, to::TimerOutput) -> LongRangeCache

Initialise the cache for the long-range backend defined in `p`.

Note that, if `pc.α === Zero()`, then long-range computations are disabled and
this returns a [`NullLongRangeCache`](@ref).
"""
function init_cache_long(pc::ParamsCommon, args...)
    if pc.α === Zero()
        NullLongRangeCache()  # disables Ewald method / long-range computations
    else
        init_cache_long_ewald(pc, args...)
    end
end

"""
    reset_fields!(cache::LongRangeCache)

Reset cache fields if required by the chosen backend.
"""
reset_fields!(c::LongRangeCache) = c

"""
    set_num_points!(cache::LongRangeCache, Np::Integer)

Set the total number of non-uniform points that the cache must hold.

Depending on the chosen [`LongRangeBackend`](@ref), this may reallocate space
to make all points fit in the cache. It will also reset to zero the
contributions of previously-added charges.

Must be called before [`add_point!`](@ref) and [`add_pointcharge!`](@ref).
"""
function set_num_points! end

"""
    add_pointcharge!(cache::LongRangeCache, X::Vec3, Q::Vec3, i::Int)

Add a vector charge at a single non-uniform location.

This is used for type-1 NUFFTs, to transform from non-uniform data in physical
space to uniform data in Fourier space.

The total number of locations must be first set via [`set_num_points!`](@ref).
"""
function add_pointcharge! end

"""
    add_point!(cache::LongRangeCache, X::Vec3, i::Int)

Add an interpolation point for type-2 NUFFT.

This is used for type-2 NUFFTs, to transform from uniform data in Fourier space
to non-uniform data in physical space.

The total number of locations must be first set via [`set_num_points!`](@ref).
"""
function add_point! end

"""
    transform_to_fourier!(cache::LongRangeCache)

Transform stored non-uniform data to Fourier space.

This usually corresponds to a type-1 NUFFT.

Non-uniform data must be first added via [`add_pointcharge!`](@ref).
"""
function transform_to_fourier! end

"""
    to_filtered_velocity!(cache::LongRangeCache)

Convert Fourier-transformed vorticity field to coarse-grained velocity field in
Fourier space.

This function should also scale the magnitude of the velocity field according
to the chosen vortex circulation `Γ` and to the cell unit dimensions `Ls`.
"""
function to_filtered_velocity!(c::LongRangeCache)
    (; uhat, ewald_op, wavenumbers,) = c
    @assert size(uhat) === size(ewald_op)
    @inbounds for I ∈ CartesianIndices(ewald_op)
        op = ewald_op[I]
        op_times_k⃗ = Vec3(map((k, i) -> @inbounds(op * k[i]), wavenumbers, Tuple(I)))
        uhat[I] = op_times_k⃗ × (im * uhat[I])
    end
    c
end

"""
    interpolate_to_physical!(cache::LongRangeCache)

Perform type-2 NUFFT to interpolate values in `cache.uhat` to non-uniform
points in physical space.

Results are written to `cache.charges`.
"""
function interpolate_to_physical! end

function _count_charges(quad::AbstractQuadrature, fs::AbstractVector{<:ClosedFilament})
    Nq = length(quad)        # number of evaluation points per filament segment
    Np = sum(f -> length(segments(f)), fs)  # total number of segments among all filaments (assumes closed filaments!!)
    Np * Nq
end

function init_ewald_fourier_operator!(
        u::AbstractArray{T, 3} where {T}, ks, Γ::Real, α::Real, Ls::NTuple{3},
    )
    prefactor = Γ / prod(Ls)
    β = -1 / (4 * α^2)
    for I ∈ CartesianIndices(u)
        k⃗ = map(getindex, ks, Tuple(I))
        k² = sum(abs2, k⃗)
        # Operator converting vorticity to coarse-grained streamfunction
        y = prefactor * exp(β * k²) / k²
        u[I] = ifelse(iszero(k²), zero(y), y)
    end
    u
end

function init_ewald_fourier_operator(::Type{T}, ks, args...) where {T <: Real}
    dims = map(length, ks)
    u = Array{T}(undef, dims)
    init_ewald_fourier_operator!(u, ks, args...)
end

function rescale_coordinates!(c::LongRangeCache)
    L = expected_period(backend(c))
    _rescale_coordinates!(c, L)
end

# Backend doesn't define `expected_period`, so no rescaling is needed.
_rescale_coordinates!(::LongRangeCache, ::Nothing) = nothing

function _rescale_coordinates!(c::LongRangeCache, L_expected::Real)
    (; Ls,) = c.params.common
    for (xs, L) ∈ zip(StructArrays.components(c.points), Ls)
        _rescale_coordinates!(xs, L, L_expected)
    end
    nothing
end

function _rescale_coordinates!(xs::AbstractVector, L::Real, L_expected)
    if L != L_expected
        xs .*= L_expected / L
    end
    nothing
end

# Note: This function must be called **after** `rescale_coordinates!`.
function fold_coordinates!(c::LongRangeCache)
    lims = folding_limits(backend(c))
    L = expected_period(backend(c))
    _fold_coordinates!(c, lims, L)
end

# Backend doesn't define `folding_limits`, so no rescaling is needed.
_fold_coordinates!(::LongRangeCache, ::Nothing, ::Any) = nothing

function _fold_coordinates!(c::LongRangeCache, lims::NTuple{2}, L::Real)
    for xs ∈ StructArrays.components(c.points)
        @inbounds for (i, x) ∈ pairs(xs)
            xs[i] = _fold_coordinate(x, lims, L)
        end
    end
    nothing
end

# We assume that L = b - a.
@inline function _fold_coordinate(x::Real, (a, b), L::Real)
    while x ≥ b
        x -= L
    end
    while x < a
        x += L
    end
    x
end

"""
    long_range_velocity_fourier!(cache::LongRangeCache, fs::AbstractVector{<:AbstractFilament})

Compute long-range velocity field induced by a set of vortex filaments.

Note that the vorticity vector averaged along the filaments must be zero to
preserve the periodicity of the velocity field (→ zero total circulation).
This function implicitly sets the average vorticity to zero.

The resulting velocity field, coarse-grained at a scale given by the Ewald
splitting parameter ``α``, is written in Fourier space to `cache.uhat`, which
is also returned by this function for convenience.
"""
function long_range_velocity_fourier!(cache::LongRangeCache, fs::VectorOfFilaments)
    (; to, params, uhat, state,) = cache
    quad = quadrature_rule(params)
    Ncharges = _count_charges(quad, fs)
    reset_fields!(cache)
    set_num_points!(cache, Ncharges)
    ζs, ws = quadrature(quad)
    n = 0
    @timeit to "add point charges" @inbounds for f ∈ fs
        ts = knots(f)
        for i ∈ eachindex(segments(f))
            Δt = ts[i + 1] - ts[i]
            for (ζ, w) ∈ zip(ζs, ws)
                X = f(i, ζ)
                Ẋ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
                # Note: the vortex circulation Γ is included in the Ewald operator and
                # doesn't need to be included here.
                q = w * Δt
                add_pointcharge!(cache, X, q * Ẋ, n += 1)
            end
        end
    end
    @assert n == Ncharges
    @timeit to "rescale coordinates" rescale_coordinates!(cache)
    @timeit to "fold coordinates" fold_coordinates!(cache)
    @timeit to "transform to Fourier" transform_to_fourier!(cache)
    state.quantity = :vorticity
    state.smoothed = false
    @timeit to "compute velocity" to_filtered_velocity!(cache)
    state.quantity = :velocity
    state.smoothed = true
    uhat
end

"""
    long_range_velocity_physical!(cache::LongRangeCache, fs::AbstractVector{<:AbstractFilament})

Interpolate long-range velocity at the location of filament nodes.

The `cache` must contain a velocity field in Fourier space. To do this, one
should first call [`long_range_velocity_fourier!`](@ref).

Velocities are written to `cache.charges`, which is also returned for convenience.
"""
function long_range_velocity_physical! end

function long_range_velocity_physical!(
        cache::LongRangeCache,
        fs::VectorOfFilaments,
    )
    (; to, charges,) = cache
    Npoints = sum(length, fs)
    set_num_points!(cache, Npoints)
    n = 0
    for f ∈ fs, X ∈ f
        add_point!(cache, X, n += 1)
    end
    @assert n == Npoints
    @timeit to "rescale coordinates" rescale_coordinates!(cache)
    @timeit to "fold coordinates" fold_coordinates!(cache)
    @timeit to "interpolate" interpolate_to_physical!(cache)
    charges
end

"""
    add_long_range_velocity!(
        vs::AbstractVector{<:AbstractVector{<:Vec3}},
        cache::LongRangeCache,
        fs::AbstractVector{<:AbstractFilament},
    )

Compute long-range part of the velocity from vortex filament locations.

The results are added to `vs`. The output array `vs` should be a vector of vectors, where
each inner vector corresponds to the nodes of a single vortex filament.
"""
function add_long_range_velocity!(
        vs::AbstractVector, cache::LongRangeCache, fs::VectorOfFilaments,
    )
    (; to,) = cache
    @timeit to "compute v̂" long_range_velocity_fourier!(cache, fs)
    @timeit to "interpolate" long_range_velocity_physical!(cache, fs)
    @timeit to "copy output" copy_long_range_output!(vs, cache)
    vs
end

function copy_long_range_output!(
        vs::AbstractVector{<:VectorOfVelocities}, cache::LongRangeCache,
    )
    (; charges,) = cache
    nout = sum(length, vs)
    nout == length(charges) || throw(DimensionMismatch("wrong length of output vector `vs`"))
    n = 0
    @inbounds for v ∈ vs, i ∈ eachindex(v)
        q = charges[n += 1]
        v[i] = v[i] + real(q)
    end
    vs
end

include("exact_sum.jl")
include("finufft.jl")
