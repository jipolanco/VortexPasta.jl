"""
    BiotSavart

Module for estimation of Biot–Savart integrals along vortex filaments using
fast Ewald splitting.
"""
module BiotSavart

export
    ParamsBiotSavart,
    GaussLegendre, NoQuadrature,
    Zero, Infinity, ∞,
    Velocity, Streamfunction,
    init_cache,
    periods,
    velocity_on_nodes!,
    compute_on_nodes!,
    reset_timer!  # from TimerOutputs

using ..BasicTypes:
    Vec3, Derivative, Zero, Infinity, ∞

using ..Quadratures:
    Quadratures, quadrature, NoQuadrature, GaussLegendre, AbstractQuadrature

using ..Filaments:
    Filaments, AbstractFilament, ClosedFilament, CurvatureBinormal,
    knots, nodes, segments, integrate

using Static: StaticBool, False, dynamic
using StructArrays: StructArrays, StructVector, StructArray
using TimerOutputs: TimerOutput, @timeit, reset_timer!

abstract type OutputField end
struct Streamfunction <: OutputField end
struct Velocity <: OutputField end

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const VectorOfVec = AbstractVector{<:Vec3}
const VectorOfPositions = VectorOfVec
const VectorOfVelocities = VectorOfVec
const AllFilamentVelocities = AbstractVector{<:VectorOfVelocities}

include("types_shortrange.jl")
include("types_longrange.jl")
include("params.jl")
include("pointdata.jl")
include("cache.jl")

include("shortrange/shortrange.jl")
include("longrange/longrange.jl")

"""
    velocity_on_nodes!(
        vs::AbstractVector{<:AbstractVector{<:Vec3}},
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament},
    ) -> vs

Compute velocity induced by vortex filaments on filament nodes.

Velocities induced by vortex filaments `fs` are written to `vs`.

This is the same as calling [`compute_on_nodes!`](@ref) when only the velocity is needed.

Usually, `fs` is a vector containing all the vortex filaments in the system.
In that case, `vs` must be a vector of vectors, which will contain the velocities of
all filament nodes. The length of `vs[i]` must be equal to the number of nodes
in the filament `fs[i]`.

The vector of velocities where the output will be written may be initialised using one of the
following lines (all are exactly equivalent):

```julia
vs = map(similar ∘ nodes, fs)
vs = [similar(nodes(f)) for f ∈ fs]
vs = similar.(nodes.(fs))
```

which initialise a velocity vector for each node of each filament (see also
[`nodes`](@ref)).
"""
function velocity_on_nodes! end

function _reset_vectors!(vs)
    for v ∈ vs
        fill!(v, zero(eltype(v)))
    end
    vs
end

function velocity_on_nodes!(
        vs::AbstractVector{<:VectorOfVec},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
        kws...,
    )
    fields = (; velocity = vs,)
    compute_on_nodes!(fields, cache, fs; kws...)
    vs
end

# Extracts and resets output fields (velocity and/or streamfunction).
function _setup_fields!(fields::NamedTuple, fs)
    N = length(fields)
    @assert N > 0
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    vecs = filter(!isnothing, (vs, ψs))
    nfields = length(vecs)
    nfields == N || throw(ArgumentError(lazy"some of these fields were not recognised: $Names"))
    for us ∈ vecs
        eachindex(us) == eachindex(fs) || throw(DimensionMismatch("wrong dimensions of vector"))
        _reset_vectors!(us)
    end
    (; vs, ψs,)
end

# Returns a tuple (Streamfunction() => ψs, Velocity() => vs) if both fields are available.
# Otherwise, it returns a subset of this which only includes the available fields.
function _fields_to_pairs(fields::NamedTuple{Names}) where {Names}
    N = length(fields)
    @assert N > 0
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    ps = (
        _make_field_pair(Velocity(), vs)...,
        _make_field_pair(Streamfunction(), ψs)...,
    )
    length(ps) == N || throw(ArgumentError(lazy"some of these fields were not recognised: $Names"))
    ps
end

_make_field_pair(key, ::Nothing) = ()
_make_field_pair(key, vs) = (key => vs,)

"""
    compute_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament};
        LIA = Val(true),
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}

Compute velocity and/or streamfunction on filament nodes.

The first argument contains one or more output fields to compute. It is usually of length 1
or 2, and can contain fields named `velocity` and `streamfunction`.

For example, to compute both velocity and streamfunction on the nodes of filaments `fs`:

```julia
# Initialise fields to compute (vectors of vectors)
vs = map(similar ∘ nodes, fs)  # one velocity vector per filament node
ψs = map(similar, vs)

# The first argument to `compute_on_nodes!` must have the following form.
# One can also choose to pass just one of the two fields.
fields = (;
    velocity = vs,
    streamfunction = ψs,
)

cache = BiotSavart.init_cache(...)
compute_on_nodes!(fields, cache, fs)
```

## Disabling LIA / computing *only* LIA

One may disable computation of the locally-induced velocity and streamfunction (LIA term)
by passing `LIA = Val(false)`. Conversely, one can pass `LIA = Val(:only)` to compute *only*
the LIA term. This can be useful for splitting the induced filament
velocities/streamfunctions onto local and non-local parts.
"""
function compute_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
        LIA = Val(true),
        longrange = Val(true),
        shortrange = Val(true),
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    if LIA === Val(:only)
        return _compute_LIA_on_nodes!(fields, cache, fs)
    end

    (; to, params, pointdata,) = cache
    (; quad,) = params
    (; vs, ψs,) = _setup_fields!(fields, fs)

    # This is used by both short-range and long-range computations.
    # Note that we need to compute the short-range first, because the long-range
    # computations then modify `pointdata`.
    @timeit to "Add point charges" add_point_charges!(pointdata, fs, quad)

    if shortrange === Val(true)
        @timeit to "Short-range component" begin
            @timeit to "Set point charges" process_point_charges!(cache.shortrange, pointdata)  # useful in particular for cell lists
            @timeit to "Compute Biot–Savart" for i ∈ eachindex(fs)
                fields_i = map(us -> us[i], fields)  # velocity/streamfunction of i-th filament
                add_short_range_fields!(fields_i, cache.shortrange, fs[i]; LIA)
            end
        end
    end

    if cache.longrange !== NullLongRangeCache() && longrange === Val(true)
        @timeit to "Long-range component" begin
            @timeit to "Vorticity to Fourier" compute_vorticity_fourier!(cache.longrange)  # uses `pointdata`
            set_interpolation_points!(cache.longrange, fs)
            if ψs !== nothing
                @timeit to "Streamfunction" begin
                    to_smoothed_streamfunction!(cache.longrange)
                    interpolate_to_physical!(cache.longrange)
                    add_long_range_output!(ψs, cache.longrange)
                end
            end
            if vs !== nothing
                # Velocity must be computed after streamfunction if both are enabled.
                @timeit to "Velocity" begin
                    to_smoothed_velocity!(cache.longrange)
                    interpolate_to_physical!(cache.longrange)
                    add_long_range_output!(vs, cache.longrange)
                end
            end
        end
    end

    fields
end

function _compute_LIA_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    (; to,) = cache
    (; params,) = cache.shortrange
    # Note: we must use the same quadrature as used when computing the globally induced terms
    (; quad, regularise_binormal,) = params
    (; Γ, a, Δ,) = params.common
    ps = _fields_to_pairs(fields)  # e.g. (Velocity() => vs, Streamfunction() => ψs)
    prefactor = Γ / (4π)
    @timeit to "LIA term (only)" begin
        @inbounds for (n, f) ∈ pairs(fs), i ∈ eachindex(f)
            for (quantity, values) ∈ ps
                # Here `quantity` is either Velocity() or Streamfunction()
                values[n][i] = BiotSavart.local_self_induced(
                    quantity, f, i, prefactor;
                    a, Δ, quad, regularise_binormal,
                )
            end
        end
    end
    fields
end

end
