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

# This may be overloaded by different backends (but it's not necessary).
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

abstract type EwaldComponent end
struct ShortRange <: EwaldComponent end
struct LongRange <: EwaldComponent end

ewald_screening_function(::Velocity, ::ShortRange, αr::Real) = erfc(αr) + 2αr / sqrt(π) * exp(-αr^2)
ewald_screening_function(::Velocity, ::ShortRange,   ::Zero) = 1

ewald_screening_function(::Velocity, ::LongRange, αr::Real) = erf(αr) - 2αr / sqrt(π) * exp(-αr^2)
ewald_screening_function(::Velocity, ::LongRange,   ::Zero) = Zero()

ewald_screening_function(::Streamfunction, ::ShortRange, αr::Real) = erfc(αr)
ewald_screening_function(::Streamfunction, ::ShortRange,   ::Zero) = 1

ewald_screening_function(::Streamfunction, ::LongRange, αr::Real) = erf(αr)
ewald_screening_function(::Streamfunction, ::LongRange,   ::Zero) = Zero()

biot_savart_integrand(::Velocity, s⃗′, r⃗, r) = (s⃗′ × r⃗) / r^3
biot_savart_integrand(::Streamfunction, s⃗′, r⃗, r) = s⃗′ / r

# Compute Biot-Savart integral over a single filament segment.
# Note: this doesn't include the prefactor Γ/4π.
function integrate_biot_savart(
        quantity::OutputField,
        component::EwaldComponent,
        seg::Segment,
        x⃗::Vec3,
        params::ParamsShortRange;
        Lhs = map(L -> L / 2, params.common.Ls),  # this allows to precompute Ls / 2
    )
    (; common, quad, rcut,) = params
    (; Ls, α,) = common
    (; f, i,) = seg
    rc² = rcut * rcut
    integrate(seg, quad) do ζ
        s⃗ = f(i, ζ)
        ∂ₜs⃗ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
        r⃗ = deperiodise_separation(x⃗ - s⃗, Ls, Lhs)
        r² = sum(abs2, r⃗)
        # Note: the short-range backend may allow some pair interactions beyond the cutoff.
        # This is in particular the case for cell lists.
        # We want to avoid extra computations in that case.
        if component === ShortRange() && r² > rc²
            zero(r⃗)
        else
            r = sqrt(r²)
            g = ewald_screening_function(quantity, component, α * r)
            w⃗ = biot_savart_integrand(quantity, ∂ₜs⃗, r⃗, r)
            oftype(r⃗, g * w⃗)  # we assume w⃗ is a vector...
        end
    end
end

"""
    add_short_range_velocity!(
        vs::AbstractVector{<:Vec3}, cache::ShortRangeCache, f::AbstractFilament;
        LIA = Val(true),
    )

Compute short-range velocity induced on the nodes of filament `f`.

The velocity vector `vs` must have the same length as the number of nodes in `f`.

Setting `LIA = Val(false)` allows to disable computation of the localised induction
approximation (LIA) term. In that case, that term should be computed separately using
[`local_self_induced_velocity`](@ref).

Before calling this function, one must first set the list of filaments using
[`set_filaments!`](@ref).
"""
function add_short_range_velocity!(
        vs::AbstractVector{<:Vec3},
        cache::ShortRangeCache,
        f::AbstractFilament;
        LIA = Val(true),  # can be used to disable LIA
    )
    fields = (; velocity = vs,)
    add_short_range_fields!(fields, cache, f; LIA)
end

function add_short_range_fields!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::ShortRangeCache,
        f::AbstractFilament;
        LIA::Val{_LIA} = Val(true),  # can be used to disable LIA
    ) where {Names, N, V <: VectorOfVec, _LIA}
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    vecs = filter(!isnothing, (vs, ψs))
    nfields = length(vecs)
    @assert nfields == N > 0

    (; params,) = cache
    (; quad,) = params
    (; Γ, a, Δ,) = params.common
    prefactor = Γ / 4π

    Xs = nodes(f)
    foreach(vecs) do us
        eachindex(us) == eachindex(Xs) || throw(DimensionMismatch(
            "vector has wrong length"
        ))
    end

    segment_a = lastindex(segments(f))   # index of segment ending at point x⃗
    segment_b = firstindex(segments(f))  # index of segment starting at point x⃗

    for (i, x⃗) ∈ pairs(Xs)
        # Start with the "singular" region (i.e. the segments which include x⃗: `segment_a` and `segment_b`).
        # We first subtract the effect of the long-range estimation, and then replace it with the LIA term.
        # Removing the long-range integral is needed to have a total velocity which does not
        # depend on the (unphysical) Ewald parameter α.
        # Note that the integral with the long-range kernel is *not* singular (since it's a
        # smoothing kernel), so there's no problem with evaluating this integral.
        # Note: we use a prefactor of 1 (instead of Γ/4π), since we intend to add the prefactor later.
        if vs !== nothing
            v⃗ = -(
                integrate_biot_savart(Velocity(), LongRange(), Segment(f, segment_a), x⃗, params) +
                integrate_biot_savart(Velocity(), LongRange(), Segment(f, segment_b), x⃗, params)
            )
            if _LIA
                v⃗ = v⃗ + local_self_induced(Velocity(), f, i, one(prefactor); a, Δ, quad,)
            end
        end
        if ψs !== nothing
            ψ⃗ = -(
                integrate_biot_savart(Streamfunction(), LongRange(), Segment(f, segment_a), x⃗, params) +
                integrate_biot_savart(Streamfunction(), LongRange(), Segment(f, segment_b), x⃗, params)
            )
            if _LIA
                ψ⃗ = ψ⃗ + local_self_induced(Streamfunction(), f, i, one(prefactor); a, Δ, quad,)
            end
        end

        # Then integrate short-range effect of nearby segments.
        for seg ∈ nearby_segments(cache, x⃗)
            g, j = seg.f, seg.i

            if f === g && (j == segment_a || j == segment_b)
                # Target point x⃗ is one of the segment limits.
                # This is the singular region, and we do nothing.
                @assert x⃗ === g[i]
            else
                # Usual case: non-singular region.
                if vs !== nothing
                    v⃗ = v⃗ + integrate_biot_savart(Velocity(), ShortRange(), seg, x⃗, params)
                end
                if ψs !== nothing
                    ψ⃗ = ψ⃗ + integrate_biot_savart(Streamfunction(), ShortRange(), seg, x⃗, params)
                end
            end
        end

        segment_a = segment_b
        segment_b += 1

        if vs !== nothing
            vs[i] = vs[i] + v⃗ * prefactor
        end
        if ψs !== nothing
            ψs[i] = ψs[i] + ψ⃗ * prefactor
        end
    end

    fields
end

include("lia.jl")  # defines local_self_induced_velocity (computation of LIA term)
include("naive.jl")
include("cell_lists.jl")
