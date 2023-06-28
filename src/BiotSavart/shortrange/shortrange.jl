using LinearAlgebra: ×, norm
using SpecialFunctions: erfc, erf
using ..Filaments: deperiodise_separation, Segment, segments

"""
    ShortRangeBackend

Abstract type denoting the backend used for computing short-range interactions.

# Implemented backends

- [`CellListsBackend`](@ref): most efficient when the cutoff radius is much smaller than the
  domain size. Can only be used with periodic boundary conditions.

- [`NaiveShortRangeBackend`](@ref): usually less efficient as it needs to compute distances
  between all pairs of locations.

# Extended help

## Implementation details

A `BACKEND <: ShortRangeBackend` must implement the function:

    init_cache_short(c::ParamsCommon, p::ParamsShortRange{<:BACKEND}, fs::AbstractVector{<:AbstractFilament}, to::TimerOutput)

which should return a [`ShortRangeCache`](@ref).
"""
abstract type ShortRangeBackend end

"""
    ShortRangeCache

Abstract type describing the storage of data required to compute short-range interactions.

The [`init_cache_short`](@ref) function returns a concrete instance of a `ShortRangeCache`.

# Interface

## Fields

The following fields must be included in a cache:

- `params <: ParamsShortRange` parameters for short-range computations;

- `to :: TimerOutput` for measuring time spent on different functions.

## Functions

The following functions must be implemented by a cache:

- [`set_filaments!`](@ref),

- [`nearby_segments`](@ref).

"""
abstract type ShortRangeCache end

"""
    init_cache_short(
        pc::ParamsCommon, p::ParamsShortRange,
        fs::AbstractVector{<:AbstractFilament},
        to::TimerOutput,
    ) -> ShortRangeCache

Initialise the cache for the short-range backend defined in `p`.
"""
function init_cache_short end

"""
    set_filaments!(c::ShortRangeCache, fs::AbstractVector{<:AbstractFilament})

Store (and optionally process) the list of filaments to be used for short-range computations.

This must be called before computing any short-range quantities (e.g. using
[`add_short_range_velocity!`](@ref)).
"""
function set_filaments! end

abstract type NearbySegmentIterator{S <: Segment} end

# These are needed e.g. by collect(it::NearbySegmentIterator)
Base.IteratorSize(::Type{<:NearbySegmentIterator}) = Base.SizeUnknown()
Base.eltype(::Type{<:NearbySegmentIterator{S}}) where {S <: Segment} = S

"""
    nearby_segments(c::ShortRangeCache, x⃗::Vec3) -> NearbySegmentIterator

Return an iterator over the segments that are "close" to the location `x⃗`.

A segment is considered to be close to `x⃗` if the minimum[^mindist] distance between `x⃗`
and the segment midpoint is smaller than the cutoff distance ``r_{\\text{cut}}``.

Typical usage:

```julia
x⃗ = Vec3(0.1, 0.3, 0.2)
for segment ∈ nearby_segments(c, x⃗)
    # Get the filament `f` and the index `i` of the segment within the filament.
    # The segment is between the filament nodes `f[i]` and `f[i + 1]`.
    (; f, i,) = segment
    # Do something with the segment...
end
```

[^mindist]: When periodicity is enabled, the relevant distance for short-range interactions
            is the *minimum* distance between the two points, after considering all their
            periodic images.
"""
function nearby_segments end

struct ParamsShortRange{
        Backend <: ShortRangeBackend,
        Quadrature <: AbstractQuadrature,
        Common <: ParamsCommon,
        T <: Real,
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    rcut    :: T           # cutoff distance

    function ParamsShortRange(
            backend::ShortRangeBackend, quad::AbstractQuadrature,
            common::ParamsCommon, rcut::Real,
        )
        (; Ls,) = common
        2 * rcut ≤ min(Ls...) ||
            error(lazy"cutoff distance `rcut = $rcut` is too large. It must be less than half the cell unit size `L` in each direction: Ls = $Ls.")
        new{typeof(backend), typeof(quad), typeof(common), typeof(rcut)}(
            backend, quad, common, rcut,
        )
    end
end

function short_range_velocity end

kernel_velocity_shortrange(αr) = erfc(αr) + 2αr / sqrt(π) * exp(-αr^2)
kernel_velocity_shortrange(::Zero) = 1

kernel_velocity_longrange(αr) = erf(αr) - 2αr / sqrt(π) * exp(-αr^2)
kernel_velocity_longrange(::Zero) = Zero()

# Compute Biot-Savart integral over a single filament segment.
# Note: this doesn't include the prefactor Γ/4π.
function integrate_biot_savart(
        kernel::F,
        seg::Segment,
        x⃗::Vec3,
        params::ParamsShortRange;
        Lhs = map(L -> L / 2, params.common.Ls),  # this allows to precompute Ls / 2
    ) where {F <: Function}
    (; common, quad,) = params
    (; Ls, α,) = common
    (; f, i,) = seg
    integrate(seg, quad) do ζ
        X = f(i, ζ)
        Ẋ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
        r⃗ = deperiodise_separation(x⃗ - X, Ls, Lhs)
        r² = sum(abs2, r⃗)
        r = sqrt(r²)
        (kernel(α * r) / r^3) * (Ẋ × r⃗)
    end
end

"""
    add_short_range_velocity!(vs::AbstractVector{<:Vec3}, cache::ShortRangeCache, f::AbstractFilament)

Compute short-range velocity induced on the nodes of filament `f`.

The velocity vector `vs` must have the same length as the number of nodes in `f`.

Before calling this function, one must first set the list of filaments using
[`set_filaments!`](@ref).
"""
function add_short_range_velocity!(
        vs::AbstractVector{<:Vec3},
        cache::ShortRangeCache,
        f::AbstractFilament;
        LIA::Val{_LIA} = Val(true),  # can be used to disable LIA (for testing only)
    ) where {_LIA}
    (; params,) = cache
    (; quad,) = params
    (; Γ, a, Δ,) = params.common
    prefactor = Γ / 4π

    Xs = nodes(f)
    eachindex(vs) == eachindex(Xs) || throw(DimensionMismatch(
        "vector of velocities has wrong length"
    ))

    segment_a = lastindex(segments(f))   # index of segment ending at point x⃗
    segment_b = firstindex(segments(f))  # index of segment starting at point x⃗

    for (i, x⃗) ∈ pairs(Xs)
        # Start with the LIA term in the singular region.
        # Note: we use a prefactor of 1 (instead of Γ/4π), since we intend to add the
        # prefactor later.
        v⃗ = local_self_induced_velocity(f, i, one(prefactor); a, Δ, quad,)
        if !_LIA
            v⃗ = zero(v⃗)
        end

        for seg ∈ nearby_segments(cache, x⃗)
            g, j = seg.f, seg.i

            # Singular region.
            if f === g && (j == segment_a || j == segment_b)
                # Target point x⃗ is one of the segment limits. This is the singular region.
                # We just need to subtract the effect of the long-range estimation, which is
                # actually included in the LIA term.
                # Without this correction is needed to have a total velocity which does not depend
                # on the (unphysical) Ewald parameter α.
                # Note that the integral with the long-range kernel is *not* singular
                # (since it's a smoothing kernel), so there's no problem with
                # evaluating this integral.
                @assert x⃗ === g[i]
                v⃗ = v⃗ - integrate_biot_savart(kernel_velocity_longrange, seg, x⃗, params)
            else
                # Usual case: non-singular region.
                v⃗ = v⃗ + integrate_biot_savart(kernel_velocity_shortrange, seg, x⃗, params)
            end
        end

        segment_a = segment_b
        segment_b += 1

        vs[i] = vs[i] + v⃗ * prefactor
    end

    vs
end

include("lia.jl")  # defines local_self_induced_velocity (computation of LIA term)
include("naive.jl")
include("cell_lists.jl")
