export CellListsBackend

using ..CellLists: CellLists, PeriodicCellList

# By default use device 2 if there are ≥ 2 available devices.
# (assuming that the long-range backend, usually NonuniformFFTsBackend, uses device 1).
shortrange_default_device(backend::KA.Backend) = min(2, KA.ndevices(backend))

@doc raw"""
    CellListsBackend <: ShortRangeBackend
    CellListsBackend([ka_backend = CPU()], [nsubdiv = 1]; [device])

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

## Using a GPU

To compute pair interactions on a GPU, pass the corresponding [KernelAbstractions.jl
backend](https://juliagpu.github.io/KernelAbstractions.jl/stable/#Supported-backends) as the
first positional argument.

For example, to use a CUDA device:

    using CUDA
    backend_short = CellListsBackend(CUDABackend(); kwargs...)

On AMD GPUs the following should work:

    using AMDGPU
    backend_short = CellListsBackend(ROCBackend(); kwargs...)

If running on a machine with multiple GPU devices, one may use the `device` keyword argument
to choose the device where short-range computations will be performed. This should be a value
in `1:ndevices`. When using KernelAbstractions.jl, the number of available devices can be
obtained using `KA.ndevices(ka_backend)`. If there are 2 or more available devices, then the
default is `device = 2` (note that `device = 1` is used by default for long-range computations,
see [`NonuniformFFTsBackend`](@ref)).

## Maximum cut-off distance

The cut-off distance must safisfy the condition:

```math
r_{\text{cut}} ≤ \frac{M}{2M + 1} L
```

where ``M`` is equal to the `nsubdiv` parameter. If this is a limitation, one can use the
[`NaiveShortRangeBackend`](@ref) which has a slightly larger limit, ``r_{\text{cut}} ≤ L/2``.
"""
struct CellListsBackend{
        M,
        BackendKA <: KA.Backend,
    } <: ShortRangeBackend
    ka_backend :: BackendKA
    ka_device  :: Int

    @inline function CellListsBackend{M}(backend::KA.Backend; device::Int) where {M}
        new{M, typeof(backend)}(backend, device)
    end
end

Base.show(io::IO, b::CellListsBackend{M}) where {M} =
    print(io, "CellListsBackend{$M}($(b.ka_backend); device = $(b.ka_device))")

@inline Base.@constprop :aggressive function CellListsBackend(
        backend::KA.Backend, n::Int = 1;
        device::Int = shortrange_default_device(backend),
    )
    CellListsBackend{n}(backend; device)
end

@inline Base.@constprop :aggressive CellListsBackend(n::Int = 1; kwargs...) =
    CellListsBackend(CPU(), n; kwargs...)

KA.get_backend(b::CellListsBackend) = b.ka_backend
KA.device(b::CellListsBackend) = b.ka_device

subdivisions(::CellListsBackend{M}) where {M} = M
max_cutoff_distance(::CellListsBackend{M}, L::AbstractFloat) where {M} = CellLists.max_cutoff_distance(M, L)

struct CellListsCache{
        Common <: ShortRangeCacheCommon,
        CellList <: PeriodicCellList,
    } <: ShortRangeCache
    common :: Common
    cl     :: CellList
end

function init_cache_short(
        pc::ParamsCommon, params::ParamsShortRange{T, <:CellListsBackend},
        pointdata::PointData, to::TimerOutput,
    ) where {T}
    (; backend, rcut,) = params
    (; Ls,) = pc
    nsubdiv = Val(subdivisions(backend))
    ka_backend = KA.get_backend(backend)
    rs_cut = map(_ -> rcut, Ls)     # same cut-off distance in each direction
    KA.device!(KA.device(backend))  # make sure we activate the targeted device
    cl = PeriodicCellList(ka_backend, rs_cut, Ls, nsubdiv)
    common = ShortRangeCacheCommon(params, pointdata, to)
    CellListsCache(common, cl)
end

function process_point_charges!(c::CellListsCache)
    (; cl, pointdata) = c
    (; points, charges, segments,) = pointdata
    @assert eachindex(points) == eachindex(charges) == eachindex(segments)
    ka_backend = KA.get_backend(c)
    @assert typeof(ka_backend) === typeof(KA.get_backend(points))
    KA.device(ka_backend) == KA.device(c) || error("expected KernelAbstractions device $(KA.device(c)) to be activated")
    Base.require_one_based_indexing(points)
    CellLists.set_elements!(cl, points)
    nothing
end

@inline nearby_charges(c::CellListsCache, x⃗::Vec3) = CellLists.nearby_elements(c.cl, x⃗)  # iterator which returns integer indices (in 1:Np)

# Note: the @inline makes a huge difference here (on Julia 1.12.1)
@inline function foreach_charge(f::F, c::CellListsCache, x⃗::Vec3; kws...) where {F <: Function}
    CellLists.foreach_source(f, c.cl, x⃗; kws...)
end

@inline function foreach_pair(f::F, c::CellListsCache; kws...) where {F <: Function}
    CellLists.foreach_pair(f, c.cl, c.pointdata.nodes; kws...)
end
