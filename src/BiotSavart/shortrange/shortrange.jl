"""
    ShortRangeBackend

Abstract type denoting the backend used for computing short-range interactions.
"""
abstract type ShortRangeBackend end

# TODO complete
"""
    ShortRangeCache

Abstract type describing the storage of data required to compute short-range interactions.
"""
abstract type ShortRangeCache end

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
