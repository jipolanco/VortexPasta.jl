using LinearAlgebra: ×
using SpecialFunctions: erfc

"""
    ShortRangeBackend

Abstract type denoting the backend used for computing short-range interactions.

# Interface

The following functions must be implemented by a `BACKEND <: ShortRangeBackend`:

- `init_cache_short(c::ParamsCommon, p::ParamsShortRange{<:BACKEND}) -> ShortRangeCache`

"""
abstract type ShortRangeBackend end

"""
    ShortRangeCache

Abstract type describing the storage of data required to compute short-range interactions.

The [`init_cache_short`](@ref) function returns a concrete instance of a `ShortRangeCache`.

# Interface

## Fields

The following fields must be included in a cache:

- `params <: ParamsShortRange` parameters for short-range computations.

"""
abstract type ShortRangeCache end

"""
    init_cache_short(pc::ParamsCommon, p::ParamsShortRange) -> ShortRangeCache

Initialise the cache for the short-range backend defined in `p`.
"""
function init_cache_short end

struct ParamsShortRange{
        Backend <: ShortRangeBackend,
        Quadrature <: AbstractQuadrature,
        Common <: ParamsCommon,
        T <: AbstractFloat,
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
        2 * rcut < min(Ls...) ||
        error(lazy"cutoff distance `rcut = $rcut` is too large. It must be less than half the cell unit size `L` in each direction: Ls = $Ls.")
        new{typeof(backend), typeof(quad), typeof(common), typeof(rcut)}(
            backend, quad, common, rcut,
        )
    end
end

# Return the shortest separation r′ given the separation r = a - b between two points given a period L.
# In the periodic line, the shortest separation satisfies |r′| ≤ L / 2.
@inline function deperiodise_separation(r::Real, L::Real, Lhalf::Real = L / 2)
    while r > Lhalf
        r -= L
    end
    while r < -Lhalf
        r += L
    end
    # @assert abs(r) ≤ Lhalf
    r
end

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

include("naive.jl")
