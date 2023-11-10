export CellListsBackend

using ..CellLists: CellLists, PeriodicCellList

"""
    CellListsBackend <: ShortRangeBackend
    CellListsBackend(nsubdiv::Int = 1)

Compute short-range interactions using the cell lists algorithm.

This backend can be significantly faster than the [`NaiveShortRangeBackend`](@ref) when the
cutoff radius `rcut` is much smaller than the domain period `L` (roughly when `rcut ≲ L / 10`).

Optionally, one can choose to subdivide each cell (of size `≈ rcut`) onto `nsubdiv`
subcells. In practice, a value of `2` or `3` can significantly improve performance compared
to no subdivision (`1`).

This backend does not support non-periodic domains.

See [`PeriodicCellList`](@ref) and [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for
more details.
"""
struct CellListsBackend{M} <: ShortRangeBackend end
CellListsBackend(n::Int = 1) = CellListsBackend{n}()

subdivisions(::CellListsBackend{M}) where {M} = M

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
        ::PointData{T}, to::TimerOutput,
    ) where {T}
    (; backend, rcut,) = params
    (; Ls,) = pc
    nsubdiv = Val(subdivisions(backend))
    @assert T <: AbstractFloat
    # Increase cut-off radius along each direction so that it exactly divides the domain size.
    rs_cut = map(Ls) do L
        L / floor(L / rcut)
    end
    @assert all(≥(rcut), rs_cut)
    Charge = Pair{Vec3{T}, Vec3{T}}       # a charge is a pair s⃗ => qs⃗′
    to_coordinate(charge) = charge.first  # retrieves the coordinate s⃗ associated to a charge
    cl = PeriodicCellList(Charge, rs_cut, Ls, nsubdiv; to_coordinate)
    CellListsCache(params, cl, to)
end

function process_point_charges!(c::CellListsCache, data::PointData)
    (; cl,) = c
    (; points, charges,) = data
    empty!(cl)
    @inbounds for i ∈ eachindex(points, charges)
        s⃗ = points[i]
        qs⃗′ = real(charges[i])  # note: `charges` may contain complex numbers (but purely real) because it's needed by some NUFFT backends
        charge = s⃗ => qs⃗′
        CellLists.add_element!(cl, charge)
    end
    nothing
end

nearby_charges(c::CellListsCache, x⃗::Vec3) = CellLists.nearby_elements(c.cl, x⃗)
