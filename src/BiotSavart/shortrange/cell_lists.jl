export CellListsBackend

using StaticArrays: StaticArrays, similar_type, Size

"""
    SegmentCellList
    SegmentCellList(::Type{<:AbstractFilament}, rcut::Real, periods::NTuple{3, Real})

Construct a cell list for dealing with the interaction between filament segments.

The applied cutoff radius ``r_{\\text{cut}}`` should be much larger than the maximum segment
length ``ℓ``, or should at least account for ``ℓ``.
Basically, if one wants an actual cut-off radius ``r₀``, then the applied cutoff radius passed
to the constructor should be ``r_{\\text{cut}} = r₀ + ℓ``.
Otherwise, a small amount of interactions within ``[r₀ - ℓ, r₀]`` may be missed.

The cutoff radius `rcut` doesn't need to exactly divide the domain period `L` into equal pieces.

Note that infinite non-periodic domains (in the sense of `period = Infinity()`) are not supported.
"""
struct SegmentCellList{
        N,  # usually N = 3 (=> 3D space)
        S <: Segment,
        CutoffRadii <: NTuple{N, Real},
        Periods <: NTuple{N, Real},
    }
    segments :: Array{Vector{S}, N}  # segments[i, j, k] contains all filament segments inside cell (i, j, k)
    rs_cut   :: CutoffRadii          # cutoff radii (can be different in each direction)
    Ls       :: Periods
end

function SegmentCellList(
        ::Type{Filament},
        rs_cut::NTuple{N, Real},
        Ls::NTuple{N, Real},
    ) where {N, Filament <: AbstractFilament}
    any(L -> L === Infinity(), Ls) && throw(ArgumentError(
        "infinite non-periodic domains not currently supported by CellListsBackend"
    ))

    # Number of cells in each direction.
    # Using `floor` below means that, if `rcut` doesn't exactly divide the domain size L in
    # a given direction, then the *last* cell in that direction will be larger than `rcut`.
    ncells = map(rs_cut, Ls) do rcut, L
        floor(Int, L / rcut)
    end

    S = Segment{Filament}
    @assert isconcretetype(S)
    segs = Array{Vector{S}, 3}(undef, ncells)
    for i ∈ eachindex(segs)
        segs[i] = S[]
    end

    SegmentCellList(segs, rs_cut, Ls)
end

function Base.empty!(cl::SegmentCellList)
    for v ∈ cl.segments
        empty!(v)
    end
    cl
end

@inline function determine_cell_index(x, rcut, L, N)
    while x < 0
        x += L
    end
    while x ≥ L
        x -= L
    end
    clamp(1 + floor(Int, x / rcut), 1, N)  # make sure the index is in 1:N
end

function add_segment!(cl::SegmentCellList{N, S}, seg::S) where {N, S <: Segment}
    (; segments, rs_cut, Ls,) = cl
    x⃗ = Filaments.midpoint(seg)
    inds = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(segments))
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

This backend can be significantly faster than the [`NaiveShortRangeBackend`](@ref) when the
cutoff radius `rcut` is much smaller than the domain period `L` (roughly when `rcut ≲ L / 10`).

Future improvements may further increase the performance of this backend.

This backend does not support non-periodic domains.

See [Wikipedia](https://en.wikipedia.org/wiki/Cell_lists) for details.
"""
struct CellListsBackend <: ShortRangeBackend end

struct CellListsCache{
        CellList <: SegmentCellList,
        Params <: ParamsShortRange,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    cl     :: CellList
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
    # Increase cut-off radius along each direction so that it exactly divides the domain size.
    rs_cut = map(Ls) do L
        L / floor(L / rcut)
    end
    @assert all(≥(rcut), rs_cut)
    cl = SegmentCellList(eltype(fs), rs_cut, Ls)
    CellListsCache(cl, params, to)
end

function set_filaments!(c::CellListsCache, fs)
    assign_cells!(c.cl, fs)
    c
end

function nearby_segments(c::CellListsCache, x⃗::Vec3)
    (; params, cl,) = c
    (; segments, rs_cut,) = cl
    (; common,) = params
    (; Ls,) = common
    inds_central = map(determine_cell_index, Tuple(x⃗), rs_cut, Ls, size(segments))
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    offsets = (-1, 0, 1)  # index offset along each dimension
    inds = map(Tuple(I₀), size(segments)) do i, N
        map(offsets) do δi
            j = i + δi
            if j ≤ 0
                j += N  # periodic wrapping
            elseif j > N
                j -= N
            end
            j  # absolute index of cell along one dimension
        end
    end
    Is = Iterators.product(inds...)  # iterate over the 3³ == 27 cells around I₀
    SA = similar_type(Array{eltype(Is)}, Size(size(Is)))  # note: the size of `Is` is known by the compiler
    Is_array = StaticArrays.sacollect(SA, Is)
    cells = map(I -> segments[I...], Is_array)  # each element is a vector of segments
    Iterators.flatten(cells)  # iterator over all segments
end
