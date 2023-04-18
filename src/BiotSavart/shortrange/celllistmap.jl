using CellListMap.PeriodicSystems:
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
