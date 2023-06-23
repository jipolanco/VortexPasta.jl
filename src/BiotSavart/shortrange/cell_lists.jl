export CellListsBackend

using ..Filaments: Segment, segments

"""
    SegmentCellList
    SegmentCellList(f::AbstractFilament, rcut::Real, periods::NTuple{3, Real})
    SegmentCellList(fs::AbstractVector{<:AbstractFilament}, rcut::Real, periods::NTuple{3, Real})

Construct a cell list for dealing with the interaction between filament segments.

In the first case, the filament `f` is only used to determine the type of segment to be
stored, and the cell list is initially empty.

In the second case, the cell list is filled with all the segments in all filaments `fs`.

The applied cutoff radius ``r_{\\text{cut}}`` should be much larger than the maximum segment
length ``ℓ``, or should at least account for ``ℓ``.
Basically, if one wants an actual cut-off radius ``r₀``, then the applied cutoff radius passed
to the constructor should be ``r_{\\text{cut}} = r₀ + ℓ``.
Otherwise, a small amount of interactions within ``[r₀ - ℓ, r₀]`` may be missed.

The cutoff radius `rcut` doesn't need to exactly divide the domain period `L` into equal pieces.

Note that infinite non-periodic domains (in the sense of `period = Infinity()`) are not supported.
"""
struct SegmentCellList{S <: Segment, T <: AbstractFloat, Periods <: NTuple{3, Real}}
    segments :: Array{Vector{S}, 3}  # segments[i, j, k] contains all filament segments inside cell (i, j, k)
    rcut     :: T                    # cutoff radius
    Ls       :: Periods
end

function SegmentCellList(f::AbstractFilament, rcut, Ls)
    any(L -> L === Infinity(), Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by CellListsBackend"
    ))

    # Number of cells in each direction (not including periodic ghost cells)
    ncells = map(Ls) do L
        ceil(Int, L / rcut)
    end

    S = eltype(segments(f))
    segs = Array{Vector{S}, 3}(undef, ncells)
    for i ∈ eachindex(segs)
        segs[i] = S[]
    end

    SegmentCellList(segs, rcut, Ls)
end

function Base.empty!(cl::SegmentCellList)
    for v ∈ cl.segments
        empty!(v)
    end
    cl
end

function SegmentCellList(fs::AbstractVector{<:AbstractFilament}, rcut, Ls)
    cl = SegmentCellList(first(fs), rcut, Ls)
    assign_cells!(cl, fs)
    cl
end

@inline function determine_cell_index(x, rcut, L)
    while x ≥ L
        x -= L
    end
    while x < 0
        x += L
    end
    1 + floor(Int, x / rcut)
end

function add_segment!(cl::SegmentCellList{S}, seg::S) where {S <: Segment}
    (; segments, rcut, Ls,) = cl
    x⃗ = Filaments.midpoint(seg)
    inds = map(Tuple(x⃗), Ls) do x, L
        @inline
        determine_cell_index(x, rcut, L)
    end
    I = CartesianIndex(inds)
    @inbounds push!(segments[I], seg)
    cl
end

function assign_cells!(cl::SegmentCellList, f::AbstractFilament)
    for s ∈ segments(f)
        add_segment!(cl, s)
    end
    cl
end

function assign_cells!(cl::SegmentCellList, fs::AbstractVector{<:AbstractFilament})
    empty!(cl)
    for f ∈ fs
        assign_cells!(cl, f)
    end
    cl
end

# ================================================================================ #

"""
    CellListsBackend <: ShortRangeBackend

Compute short-range interactions using the cell lists algorithm.

This backend does not support non-periodic domains.

See [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for details.
"""
struct CellListsBackend <: ShortRangeBackend end

struct CellListsCache{
        Params <: ParamsShortRange,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    cl     :: SegmentCellList
    params :: Params
    to     :: Timer
end

function init_cache_short(
        pc::ParamsCommon, params::ParamsShortRange{<:CellListsBackend},
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    )
    (; rcut,) = params
    (; Ls,) = pc
    cl = SegmentCellList(fs, rcut, Ls)
    CellListsCache(cl, params, to)
end
