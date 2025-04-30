export CellListsBackend

using ..CellLists: CellLists, PeriodicCellList

@doc raw"""
    CellListsBackend <: ShortRangeBackend
    CellListsBackend(nsubdiv::Int = 1)

Compute short-range interactions using the cell lists algorithm.

This backend can be significantly faster than the [`NaiveShortRangeBackend`](@ref) when the
cut-off distance `r_cut` is much smaller than the domain period `L` (roughly when `r_cut ≲ L / 10`).

Optionally, one can choose to subdivide each cell (of size `≈ r_cut`) onto `nsubdiv`
subcells. In practice, a value of `2` or `3` can significantly improve performance compared
to no subdivision (`1`).

Note that, with this backend, the cut-off distance must satisfy `r_cut ≤ M / (2M + 1) * L`
where `M = nsubdiv`.

This backend does not support non-periodic domains.

See [`PeriodicCellList`](@ref) and [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for
more details.

## Maximum cut-off distance

The cut-off distance must safisfy the condition:

```math
r_{\text{cut}} ≤ \frac{M}{2M + 1} L
```

where ``M`` is equal to the `nsubdiv` parameter. If this is a limitation, one can use the
[`NaiveShortRangeBackend`](@ref) which has a slightly larger limit, ``r_{\text{cut}} ≤ L/2``.
"""
struct CellListsBackend{M} <: ShortRangeBackend end
CellListsBackend(n::Int = 1) = CellListsBackend{n}()

subdivisions(::CellListsBackend{M}) where {M} = M
max_cutoff_distance(::CellListsBackend{M}, L::AbstractFloat) where {M} = oftype(L, M / (2M + 1)) * L

struct CellListsCache{
        Params <: ParamsShortRange{<:Real, <:CellListsBackend},
        CellList <: PeriodicCellList,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    params :: Params
    cl     :: CellList
    to     :: Timer
end

function init_cache_short(
        pc::ParamsCommon, params::ParamsShortRange{T, <:CellListsBackend},
        pdata::PointData{T}, to::TimerOutput,
    ) where {T}
    (; backend, rcut,) = params
    (; Ls,) = pc
    nsubdiv = Val(subdivisions(backend))
    @assert T <: AbstractFloat
    rs_cut = map(_ -> rcut, Ls)  # same cut-off distance in each direction
    Seg = eltype(pdata.segments)
    @assert Seg <: Segment
    Charge = Tuple{Vec3{T}, Vec3{T}, Seg}  # a "charge" is a tuple (s⃗, qs⃗′, segment)
    @assert isconcretetype(Charge)
    to_coordinate(charge) = first(charge)  # retrieves the coordinate s⃗ associated to a charge
    cl = PeriodicCellList(Charge, rs_cut, Ls, nsubdiv; to_coordinate)
    CellListsCache(params, cl, to)
end

function process_point_charges!(c::CellListsCache, data::PointData)
    (; cl,) = c
    (; points, charges, segments,) = data
    empty!(cl)
    @inbounds for i ∈ eachindex(points, charges, segments)
        s⃗ = points[i]
        qs⃗′ = real(charges[i])  # note: `charges` may contain complex numbers (but purely real) because it's needed by some NUFFT backends
        seg = segments[i]
        charge = (s⃗, qs⃗′, seg)
        CellLists.add_element!(cl, charge)
    end
    CellLists.finalise_cells!(cl)
    nothing
end

nearby_charges(c::CellListsCache, x⃗::Vec3) = CellLists.nearby_elements(c.cl, x⃗)
