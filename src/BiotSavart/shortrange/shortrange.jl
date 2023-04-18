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
        T <: AbstractFloat,
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    rcut    :: T           # cutoff distance
end

include("celllistmap.jl")
