using LinearAlgebra: ×, norm
using SpecialFunctions: erfc, erf
using ..Filaments: deperiodise_separation

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

include("lia.jl")  # defines local_self_induced_velocity (computation of LIA term)
include("naive.jl")
