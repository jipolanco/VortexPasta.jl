export CellListSegmentFinder

using ..CellLists: CellLists, PeriodicCellList, static, StaticInt

"""
    CellListSegmentFinder <: NearbySegmentFinder
    CellListSegmentFinder(
        fs::AbstractVector{<:AbstractFilament},
        r_cut::Real,
        Ls::NTuple{3, Real};
        nsubdiv = Val(1),
    )

Initialise nearby segment finder based on cell lists algorithm.

Only supports fully periodic domains.

See the [`CellLists`](@ref) module for more details.

# Mandatory arguments

- `fs`: vector of filaments;
- `r_cut::Real`: cut-off distance;
- `Ls::NTuple{3, Real}`: domain period in each direction.
  Note that this backend doesn't support non-periodic (`Infinity`) directions.

# Optional keyword arguments

- `nsubdiv = Val(M)`: number of cell lists subdivisions.
  Here `M` must be a positive integer (`M = 1` means no subdivisions).
  See [`PeriodicCellList`](@ref) for details.
    
"""
struct CellListSegmentFinder{
        CellList <: PeriodicCellList,
        CutoffDistance <: Real,
        Periods <: Tuple{Vararg{Real}},
    } <: NearbySegmentFinder
    cl     :: CellList
    r_cut  :: CutoffDistance
    Ls     :: Periods
end

function CellListSegmentFinder(
        fs::AbstractVector{<:AbstractFilament}, r_cut, Ls;
        nsubdiv = static(1),
    )
    # Increase cut-off radius along each direction so that it exactly divides the domain size.
    rs_cut = map(Ls) do L
        L / floor(L / r_cut)
    end
    @assert all(≥(r_cut), rs_cut)
    Filament = eltype(fs)
    S = Segment{Filament}
    @assert isconcretetype(S)
    # Construct cell list of filament segments.
    # We use the `midpoint` function to associate a coordinate to each segment.
    M = to_static(nsubdiv)
    cl = PeriodicCellList(S, rs_cut, Ls, M; to_coordinate = Filaments.midpoint)
    CellListSegmentFinder(cl, r_cut, Ls)
end

to_static(::Val{M}) where {M} = static(M)
to_static(M::StaticInt) = M

function _assign_cells!(cl::PeriodicCellList, f::AbstractFilament)
    for s ∈ segments(f)
        CellLists.add_element!(cl, s)
    end
    cl
end

function _assign_cells!(cl::PeriodicCellList, fs::AbstractVector{<:AbstractFilament})
    empty!(cl)
    for f ∈ fs
        _assign_cells!(cl, f)
    end
    cl
end

function set_filaments!(c::CellListSegmentFinder, fs)
    _assign_cells!(c.cl, fs)
    c
end

nearby_segments(c::CellListSegmentFinder, x⃗::Vec3) = CellLists.nearby_elements(c.cl, x⃗)