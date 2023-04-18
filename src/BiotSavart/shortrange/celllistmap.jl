using CellListMap.PeriodicSystems:
    AbstractPeriodicSystem,
    PeriodicSystem, map_pairwise!

"""
    CellListMapBackend <: ShortRangeBackend

Compute short-range interactions using the
[CellListMap.jl](https://github.com/m3g/CellListMap.jl) package.

CellListMap.jl provides computation of pair interactions in periodic systems
using cell list algorithms.

Computations can be parallelised using threads.
"""
struct CellListMapBackend <: ShortRangeBackend end

struct CellListMapCache{
        Params <: ParamsShortRange,
        System <: AbstractPeriodicSystem,
    } <: ShortRangeCache
    params :: Params
    system :: System
end

function init_cache_short(
        common::ParamsCommon{T}, params::ParamsShortRange{<:CellListMapBackend},
    ) where {T}
    (; Ls,) = common
    (; rcut,) = params
    positions = [zero(Vec3{T})]  # CellListMap needs at least one element
    system = PeriodicSystem(
        xpositions = positions,
        unitcell = SVector(Ls),
        cutoff = rcut,
        output = similar(positions),
    )
    CellListMapCache(
        params, system,
    )
end
