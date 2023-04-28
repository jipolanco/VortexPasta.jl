using LinearAlgebra: ×, norm
using SpecialFunctions: erfc, erf

"""
    ShortRangeBackend

Abstract type denoting the backend used for computing short-range interactions.

# Interface

The following functions must be implemented by a `BACKEND <: ShortRangeBackend`:

- `init_cache_short(c::ParamsCommon, p::ParamsShortRange{<:BACKEND}, to::TimerOutput) -> ShortRangeCache`

"""
abstract type ShortRangeBackend end

"""
    ShortRangeCache

Abstract type describing the storage of data required to compute short-range interactions.

The [`init_cache_short`](@ref) function returns a concrete instance of a `ShortRangeCache`.

# Interface

## Fields

The following fields must be included in a cache:

- `params <: ParamsShortRange` parameters for short-range computations;

- `to :: TimerOutput` for measuring time spent on different functions.

"""
abstract type ShortRangeCache end

"""
    init_cache_short(pc::ParamsCommon, p::ParamsShortRange, to::TimerOutput) -> ShortRangeCache

Initialise the cache for the short-range backend defined in `p`.
"""
function init_cache_short end

struct ParamsShortRange{
        Backend <: ShortRangeBackend,
        Quadrature <: AbstractQuadrature,
        Common <: ParamsCommon,
        T <: Real,
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    rcut    :: T           # cutoff distance

    function ParamsShortRange(
            backend::ShortRangeBackend, quad::AbstractQuadrature,
            common::ParamsCommon, rcut::Real,
        )
        (; Ls,) = common
        2 * rcut ≤ min(Ls...) ||
            error(lazy"cutoff distance `rcut = $rcut` is too large. It must be less than half the cell unit size `L` in each direction: Ls = $Ls.")
        new{typeof(backend), typeof(quad), typeof(common), typeof(rcut)}(
            backend, quad, common, rcut,
        )
    end
end

# Return the shortest separation r′ given the separation r = a - b between two points given a period L.
# In the periodic line, the shortest separation satisfies |r′| ≤ L / 2.
@inline function deperiodise_separation(r::Real, L::Real, Lhalf::Real)
    while r > Lhalf
        r -= L
    end
    while r < -Lhalf
        r += L
    end
    # @assert abs(r) ≤ Lhalf
    r
end

# Infinite period -> nothing to deperiodise.
@inline deperiodise_separation(r::Real, ::Infinity, ::Infinity) = r

# We convert SVector to tuple to make sure that no heap allocations are performed.
@inline deperiodise_separation(r⃗::Vec3, args...) = oftype(r⃗, deperiodise_separation(Tuple(r⃗), args...))
@inline deperiodise_separation(rs::Tuple, args...) = map(deperiodise_separation, rs, args...) :: Tuple

"""
    short_range_velocity(cache::ShortRangeCache, x⃗::Vec3, f::AbstractFilament, [inds = eachindex(segments(f))])

Compute short-range velocity induced by filament `f` on coordinate `x⃗`.

This involves the estimation of a line integral over the filament `f`.
By default, the integration is performed over the whole filament.

One may optionally choose to integrate over a subset of the segments of the filament.
This is useful for avoiding the Biot–Savart singularity when the point `x⃗`
belongs to the filament.
"""
function short_range_velocity end

kernel_velocity_shortrange(αr) = erfc(αr) + 2αr / sqrt(π) * exp(-αr^2)
kernel_velocity_shortrange(::Zero) = 1

kernel_velocity_longrange(αr) = erf(αr) - 2αr / sqrt(π) * exp(-αr^2)
kernel_velocity_longrange(::Zero) = Zero()

function short_range_velocity(cache::ShortRangeCache, args...)
    short_range_velocity(kernel_velocity_shortrange, cache, args...)
end

"""
    add_short_range_velocity_self!(vs::AbstractVector{<:Vec3}, cache::ShortRangeCache, f::AbstractFilament)

Compute short-range self-induced velocity of a filament on its own nodes.

The result is added to existent values in the `vs` vector.

This includes:

- the LIA term (*localised induction approximation*), i.e. the local self-induced
  velocity based on the local filament curvature;

- the non-local (and short-range in the Ewald sense) velocity induced by the filament.

This does *not* include:

- the short-range velocity induced by other filaments (see
  [`add_short_range_velocity_other!`](@ref) for that);

- the long-range velocity induced by all filaments (see
  [`add_long_range_velocity!`](@ref) for that).

The length of the output vector `vs` must be equal to the number of nodes of the filament `f`.
"""
function add_short_range_velocity_self!(
        vs::VectorOfVelocities,
        cache::ShortRangeCache,
        f::AbstractFilament;
        LIA::Bool = true,  # allows disabling LIA for testing
    )
    (; to, params,) = cache
    (; quad,) = params
    (; Γ, a, Δ) = params.common

    prefactor = Γ / 4π

    length(vs) == length(f) || throw(DimensionMismatch("wrong length of output `vs`"))
    N = length(segments(f))
    inds_all = UnitRange(eachindex(segments(f)))

    # Case of the initial node (i = 1): we don't integrate over the first nor
    # the last segment (assuming periodicity / closed filament...).
    inds_left = inds_all[0:-1]
    inds_right = inds_all[2:end - 1]
    @assert inds_all === 1:N
    inds_singular = 0:1  # singularity region

    for (i, x⃗) ∈ pairs(nodes(f))
        v⃗ = zero(eltype(vs))
        @assert length(inds_left) + length(inds_right) == length(inds_all) - 2
        v⃗ = v⃗ + short_range_velocity(cache, x⃗, f, inds_left)  # this already includes the prefactor
        # We subtract part of the long-range velocity computed by the long-range backend.
        # Namely, we subtract the local self-induced velocity in the
        # singularity region, which will be replaced by the LIA approximation.
        # This correction is needed to have a total velocity which does not depend
        # on the (unphysical) Ewald parameter α.
        # Note that the integral with the long-range kernel is *not* singular
        # (since it's a smoothing kernel), so there's no problem with
        # evaluating this integral.
        v⃗ = v⃗ - short_range_velocity(kernel_velocity_longrange, cache, x⃗, f, inds_singular)
        v⃗ = v⃗ + short_range_velocity(cache, x⃗, f, inds_right)
        if LIA
            v⃗_LIA = local_self_induced_velocity(f, i, prefactor; a, Δ, quad)
            v⃗ = v⃗ + v⃗_LIA
        end
        vs[i] = vs[i] + v⃗
        # For the next point (i + 1), the singularity region corresponds to the
        # segments [i, i + 1].
        inds_left = inds_all[begin:(i - 1)]
        inds_right = inds_all[(i + 2):end]
        inds_singular = i:(i + 1)
    end

    vs
end

"""
    add_short_range_velocity_other!(
        vs::AbstractVector{<:Vec3},
        Xs::AbstractVector{<:Vec3},
        cache::ShortRangeCache,
        f::AbstractFilament,
    )

Compute short-range velocity induced by a filament `f` at locations `Xs`.

The locations `Xs` can be the nodes of a vortex filament different from `f`.
They can also be arbitrary locations in the domain.
"""
function add_short_range_velocity_other!(
        vs::VectorOfVelocities,
        Xs::VectorOfPositions,
        cache::ShortRangeCache,
        f::AbstractFilament,
    )
    eachindex(vs) == eachindex(Xs) || throw(DimensionMismatch("wrong length of output `vs`"))
    for (i, x⃗) ∈ pairs(Xs)
        v⃗ = short_range_velocity(cache, x⃗, f)  # already includes the prefactor Γ/4π
        vs[i] = vs[i] + v⃗
    end
    vs
end

"""
    local_self_induced_velocity(f::AbstractFilament, i::Int; Γ::Real, a::Real, Δ = 0.25, quad = nothing)
    local_self_induced_velocity(f::AbstractFilament, i::Int, prefactor::Real; a::Real, Δ = 0.25, quad = nothing)

Compute local self-induced velocity of filament node `f[i]`.

This corresponds to the LIA term (localised induction approximation).

# Mandatory arguments

- the vortex core size `a` as a keyword argument;

- either the vortex circulation `Γ` as a keyword argument, or the precomputed
  coefficient `Γ / 4π` as a positional argument.

# Optional arguments

- the optional parameter `Δ` sets the effect of the vorticity profile (see
  [`ParamsBiotSavart`](@ref) for details);

- a quadrature rule can be passed via `quad`, which can improve the estimation
  of the LIA term for relatively large segments.
"""
function local_self_induced_velocity end

function _local_self_induced_velocity(
        quad::Nothing, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real,
    )
    ts = knots(f)
    ℓ₋ = ts[i] - ts[i - 1]  # assume that the parametrisation roughly corresponds to the vortex arc length
    ℓ₊ = ts[i + 1] - ts[i]
    β = prefactor * (log(2 * sqrt(ℓ₋ * ℓ₊) / a) - Δ)
    Ẋ = f[i, Derivative(1)]  # ∂f/∂t at node i
    Ẍ = f[i, Derivative(2)]
    # Note that the derivatives are wrt the parameter `t` and not exactly wrt
    # the arc length `ξ`, hence the extra `norm(Ẋ)^3` factor wrt the literature.
    # Usually, `norm(Ẋ)` should be actually quite close to 1 since the
    # parametrisation is a rough approximation of the arc length.
    β / norm(Ẋ)^3 * (Ẋ × Ẍ)
end

# Alternative estimation using quadratures.
# It seems to improve accuracy (tested with vortex ring example).
function _local_self_induced_velocity(
        quad::AbstractQuadrature, f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real,
    )
    ts = knots(f)
    ℓ₋ = integrate(f, i - 1, quad) do ζ
        norm(f(i - 1, ζ, Derivative(1)))
    end
    ℓ₊ = integrate(f, i, quad) do ζ
        norm(f(i, ζ, Derivative(1)))
    end
    b⃗₋ = integrate(f, i - 1, quad) do ζ
        Ẋ = f(i - 1, ζ, Derivative(1))
        Ẍ = f(i - 1, ζ, Derivative(2))
        (Ẋ × Ẍ) ./ norm(Ẋ)^3
    end
    b⃗₊ = integrate(f, i, quad) do ζ
        Ẋ = f(i, ζ, Derivative(1))
        Ẍ = f(i, ζ, Derivative(2))
        (Ẋ × Ẍ) ./ norm(Ẋ)^3
    end
    b⃗ = (b⃗₋ + b⃗₊) ./ (ts[i + 1] - ts[i - 1])  # average on [i - 1, i + 1]
    β = prefactor * (log(2 * sqrt(ℓ₋ * ℓ₊) / a) - Δ)
    β * b⃗
end

function local_self_induced_velocity(
        f::AbstractFilament, i::Int, prefactor::Real;
        a::Real, Δ::Real = 0.25, quad = nothing,
    )
    _local_self_induced_velocity(quad, f, i, prefactor; a, Δ)
end

local_self_induced_velocity(f, i; Γ, kws...) =
    local_self_induced_velocity(f, i, Γ / 4π; kws...)

include("naive.jl")
