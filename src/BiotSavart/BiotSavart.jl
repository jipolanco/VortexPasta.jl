"""
    BiotSavart

Module for estimation of Biot–Savart integrals along vortex filaments using
fast Ewald splitting.
"""
module BiotSavart

export
    ParamsBiotSavart,
    GaussLegendre, NoQuadrature, AdaptiveTanhSinh,
    Zero, Infinity, ∞,
    Velocity, Streamfunction,
    LongRangeCache, ShortRangeCache,
    init_cache,
    has_real_to_complex,
    periods,
    velocity_on_nodes!,
    compute_on_nodes!,
    CPU,  # from KernelAbstractions
    reset_timer!  # from TimerOutputs

using ..BasicTypes:
    Vec3, Derivative, Zero, Infinity, ∞

using ..Quadratures:
    Quadratures, quadrature, NoQuadrature, GaussLegendre, AdaptiveTanhSinh,
    AbstractQuadrature, StaticSizeQuadrature, PreallocatedQuadrature

using ..Filaments:
    Filaments, AbstractFilament, ClosedFilament, Segment, CurvatureBinormal,
    knots, nodes, segments, integrate

using Adapt: Adapt, adapt

using KernelAbstractions:
    KernelAbstractions,  # importing this avoids docs failure
    KernelAbstractions as KA, @kernel, @index, @Const,
    CPU, GPU

using StableTasks: StableTasks

using Bumper: Bumper, @no_escape, @alloc
using ChunkSplitters: ChunkSplitters
using StructArrays: StructArrays, StructVector, StructArray
using TimerOutputs: TimerOutputs, TimerOutput, @timeit, reset_timer!

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
include("autotune.jl")

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
function _setup_fields!(fields::NamedTuple{Names}, fs) where {Names}
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
# If `i` is an integer, then `ψs` and `vs` only contain the data associated to a single filament `i`.
Base.@propagate_inbounds function _fields_to_pairs(fields::NamedTuple{Names}, i = nothing) where {Names}
    N = length(fields)
    @assert N > 0
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    ps = (
        _make_field_pair(Velocity(), vs, i)...,
        _make_field_pair(Streamfunction(), ψs, i)...,
    )
    length(ps) == N || throw(ArgumentError(lazy"some of these fields were not recognised: $Names"))
    ps
end

_make_field_pair(key, ::Nothing, ::Any) = ()
_make_field_pair(key, vs::AbstractVector, ::Nothing) = (key => vs,)
Base.@propagate_inbounds _make_field_pair(key, vs::AbstractVector, i::Integer) = (key => vs[i],)

@doc raw"""
    background_vorticity_correction!(
        fields::NamedTuple, fs::AbstractVector{<:AbstractFilament}, params::ParamsBiotSavart,
    )

Correct computed fields in case the mean filament vorticity is non-zero.

This ensures that results do not depend on the Ewald splitting parameter ``α`` when the
filament "charge" is non-zero in a periodic domain.

In practice, as detailed below, this function only corrects the short-range streamfunction,
as the velocity is not affected by the background vorticity, and long-range corrections are
done implicitly in Fourier space.

!!! note

    This function is automatically called by [`compute_on_nodes!`](@ref) and should never be
    called manually.
    This documentation is only included for information purposes.

# Extended help

The mean vorticity associated to the filaments is given by

```math
⟨ \bm{ω} ⟩ = -\frac{Γ}{V} ∮_{\mathcal{C}} \mathrm{d}\bm{s},
```

where $V = L^3$ is the periodic domain volume.

This quantity is always zero when all filaments are closed. It can be non-zero if there are
infinite unclosed filaments which are not compensated by oppositely-oriented infinite filaments.

When the mean filament vorticity is non-zero, one must compensate it by adding a uniform background
vorticity $\bm{ω}_{\text{back}} = -⟨ \bm{ω} ⟩$, so that the total circulation around the
volume is zero.

In practice, in the long-range component of Ewald summation, this compensation is done
implicitly by setting the zero mode $\hat{\bm{\omega}}(\bm{k} = \bm{0}) = \bm{0}$. As for
the short-range component, only the streamfunction needs to be corrected, while the velocity
is not affected by the background vorticity.

The correction to the short-range streamfunction is given by

```math
\bm{ψ}^{<}_{\text{back}}
= \frac{\bm{ω}_{\text{back}}}{4π} ∫_{\mathbb{R}^3} \frac{\operatorname{erfc}(α|\bm{r}|)}{|\bm{r}|} \, \mathrm{d}^3\bm{r}
= \frac{\bm{ω}_{\text{back}}}{4 α^2}
```

Therefore, this function simply adds this short-range correction to the streamfunction.

"""
function background_vorticity_correction!(
        fields::NamedTuple,
        fs::AbstractVector{<:AbstractFilament},
        params::ParamsBiotSavart,
    )
    # We only need to modify the streamfunction, so we do nothing if the streamfunction is
    # not present. (The velocity is not affected by the uniform background vorticity.)
    ψs_all = get(fields, :streamfunction, nothing)
    ψs_all === nothing && return fields
    domain_is_periodic(params) || return fields  # correction only makes sense in periodic domain

    (; Γ, Ls, α,) = params

    # Compute mean vorticity associated to filaments.
    # This simply corresponds to adding the end_to_end_offset's of each filament.
    ω⃗_mean = sum(Filaments.end_to_end_offset, fs) * (Γ / prod(Ls))
    iszero(ω⃗_mean) && return fields  # there's nothing to do if the mean vorticity is zero

    ω⃗_back = -ω⃗_mean  # background vorticity

    # Streamfunction associated to background vorticity
    ψ⃗_back = ω⃗_back / (4 * α^2)

    @inbounds for ψs ∈ ψs_all, i ∈ eachindex(ψs)
        ψs[i] = ψs[i] + ψ⃗_back
    end

    fields
end

"""
    compute_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament};
        LIA = Val(true),
        shortrange = true,
        longrange = true,
        callback_vorticity = identity,
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

# Extended help

## Disabling local terms / computing *only* local terms

One may disable computation of the locally-induced velocity and streamfunction (LIA term)
by passing `LIA = Val(false)`. Conversely, one can pass `LIA = Val(:only)` to compute *only*
the LIA term. This can be useful for splitting the induced filament
velocities/streamfunctions onto local and non-local parts.

## Disabling short-range or long-range interactions

It is also possible to disable the short-range or long-range component of Ewald splitting,
if one only wants to compute one of the two components.
To do this, pass either `shortrange = false` or `longrange = false`.

Note that the short-range component includes the local (LIA) term as well as the short-range
correction term for long-range interactions. Therefore, setting `shortrange = false`
disables both these terms.

## Accessing vorticity in Fourier space

The computation of long-range quantities involves estimating the Fourier coefficients of the
vorticity field associated to the vortex filaments. These coefficients are truncated to some
maximum wavenumber ``k_{\\text{max}}`` in each Cartesian direction. This information can be
useful for other things, for instance computing energy spectra.

One can use the `callback_vorticity` argument to access the vorticity in Fourier space,
before it is replaced by the coefficients of streamfunction and/or velocity.
This argument should be a function `callback_vorticity(cache)` which takes a
[`LongRangeCache`](@ref). The callback should not modify anything inside the cache, or
otherwise the streamfunction and velocity computed by this function will likely be wrong.
Of course, this callback will be ignored if long-range computations are disabled.

Note that, when the callback is called, the vorticity coefficients in `cache.common.uhat`
don't have the right physical dimensions as they have not yet been multiplied by
``Γ/V`` (where ``V`` is the volume of a unit cell). Note that ``Γ/V`` is directly available
in `cache.common.ewald_prefactor`. Besides, the vorticity coefficients at this stage have
not yet been Gaussian-filtered according to Ewald's method.

An example of how to compute the (large-scale) kinetic energy associated to the
Fourier-truncated vorticity field:

```julia
using Adapt: adapt  # useful in case FFTs are computed on the GPU

E_from_vorticity = Ref(0.0)  # "global" variable updated when calling compute_on_nodes!

function callback_vorticity(cache::LongRangeCache)
    (; wavenumbers_d, uhat_d, ewald_prefactor,) = cache.common
    # For simplicity, copy data to the CPU if it's on the GPU.
    wavenumbers = adapt(Array, wavenumbers_d)
    uhat = adapt(Array, uhat_d)
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)  # this depends on the long-range backend
    @assert with_hermitian_symmetry == wavenumbers[1][end] > 0
    γ² = ewald_prefactor^2  # = (Γ/V)^2 [prefactor not included in the vorticity]
    E = 0.0
    for I ∈ CartesianIndices(uhat)
        k⃗ = map(getindex, wavenumbers, Tuple(I))
        kx = k⃗[1]
        factor = (!with_hermitian_symmetry || kx == 0) ? 0.5 : 1.0
        k² = sum(abs2, k⃗)
        if !iszero(k²)
            ω⃗ = uhat[I]  # Fourier coefficient of the vorticity
            E += γ² * factor * sum(abs2, ω⃗) / k²
        end
    end
    E_from_vorticity[] = E  # update value of "global" variable
    nothing
end

compute_on_nodes!(...; callback_vorticity)
```

"""
function compute_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
        LIA = Val(true),
        kws...,
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    if LIA === Val(:only)
        return _compute_LIA_on_nodes!(fields, cache, fs)
    end
    device_lr = KA.get_backend(cache.longrange)  # device used for long-range part (CPU or GPU)
    _compute_on_nodes!(device_lr, fields, cache, fs; LIA, kws...)
    fields
end

# CPU/CPU version (no GPU)
function _compute_on_nodes!(
        ::CPU, fields::NamedTuple, cache, fs;
        LIA = Val(true),
        longrange = true,
        shortrange = true,
        callback_vorticity::Fvort = identity,
    ) where {Fvort}
    (; to, params, pointdata,) = cache
    (; quad,) = params
    (; vs, ψs,) = _setup_fields!(fields, fs)

    with_shortrange = shortrange
    with_longrange = longrange && cache.longrange !== NullLongRangeCache()

    # This is used by both short-range and long-range computations.
    # Note that we need to compute the short-range first, because the long-range
    # computations then modify `pointdata`.
    @timeit to "Add point charges" add_point_charges!(pointdata, fs, quad)  # on the CPU

    if with_shortrange
        @timeit to "Short-range component" begin
            @timeit to "Process point charges" process_point_charges!(cache.shortrange, pointdata)  # useful in particular for cell lists
            @timeit to "Compute Biot–Savart" add_short_range_fields!(fields, cache.shortrange, fs; LIA)
            @timeit to "Background vorticity" background_vorticity_correction!(fields, fs, params)
        end
    end

    if with_longrange
        @timeit to "Long-range component" begin
            @timeit to "Vorticity to Fourier" begin
                compute_vorticity_fourier!(cache.longrange)  # reads pointdata (points and charges)
            end
            if callback_vorticity !== identity
                @timeit to "Vorticity callback" begin
                    callback_vorticity(cache.longrange)
                end
            end
            @timeit to "Set interpolation points" begin
                set_interpolation_points!(cache.longrange, fs)  # overwrites pointdata (points)
            end
            if ψs !== nothing
                @timeit to "Streamfunction" begin
                    @timeit to "Convert to physical" begin
                        to_smoothed_field!(Streamfunction(), cache.longrange)
                        interpolate_to_physical!(cache.longrange)  # overwrites pointdata (charges)
                        add_long_range_output!(ψs, cache.longrange)
                    end
                    @timeit to "Self-interaction" remove_long_range_self_interaction!(
                        ψs, fs, Streamfunction(), params.common,
                    )
                end
            end
            if vs !== nothing
                # Velocity must be computed after streamfunction if both are enabled.
                @timeit to "Velocity" begin
                    @timeit to "Convert to physical" begin
                        to_smoothed_field!(Velocity(), cache.longrange)
                        interpolate_to_physical!(cache.longrange)  # overwrites pointdata (charges)
                        add_long_range_output!(vs, cache.longrange)
                    end
                    @timeit to "Self-interaction" remove_long_range_self_interaction!(
                        vs, fs, Velocity(), params.common,
                    )
                end
            end
        end
    end

    nothing
end

# CPU/GPU version
# We compute short-range (CPU) and long-range (GPU) asynchronously, so that both components
# work at the same time.
function _compute_on_nodes!(
        device_lr::GPU, fields::NamedTuple, cache, fs;
        LIA = Val(true),
        longrange = true,
        shortrange = true,
        callback_vorticity::Fvort = identity,
    ) where {Fvort}
    (; to, params, pointdata,) = cache
    (; pointdata_d, to_d,) = cache.longrange.common  # pointdata on the device (GPU)
    (; quad,) = params
    (; vs, ψs,) = _setup_fields!(fields, fs)

    nfields = length(fields)
    @assert length(fields) ∈ (1, 2)  # maximum 2 fields: streamfunction + velocity
    @assert (vs === nothing) + (ψs === nothing) + nfields == 2

    with_shortrange = shortrange
    with_longrange = longrange && cache.longrange !== NullLongRangeCache()

    # This is used by short and long-range computations.
    @timeit to "Add point charges" add_point_charges!(pointdata, fs, quad)  # done on the CPU

    # Allocate temporary arrays on the GPU for interpolation outputs (manually deallocated later).
    if with_longrange
        noutputs = sum(length, fs)  # total number of interpolation points
        outputs_lr = map(_ -> similar(pointdata_d.charges, noutputs), values(fields))
    end

    # Compute long-range part asynchronously on the GPU.
    task_lr = if with_longrange
        StableTasks.@spawn begin
            @timeit to_d "Long-range component (GPU)" begin
                # Copy point data to the GPU (pointdata_d).
                @assert pointdata !== pointdata_d  # they are different objects
                @assert device_lr isa PseudoGPU || typeof(pointdata) !== typeof(pointdata_d)
                @timeit to_d "Copy point charges (host → device)" begin
                    copy!(pointdata_d, pointdata)  # H2D copy
                    KA.synchronize(device_lr)
                end
                @timeit to_d "Vorticity to Fourier" begin
                    compute_vorticity_fourier!(cache.longrange)  # reads pointdata_d (points and charges)
                end
                if callback_vorticity !== identity
                    @timeit to_d "Vorticity callback" begin
                        callback_vorticity(cache.longrange)
                    end
                end
                @timeit to_d "Set interpolation points" begin
                    set_interpolation_points!(cache.longrange, fs)  # overwrites pointdata_d (points)
                end
                # Interpolate streamfunction and/or velocity.
                local ifield = 0
                if ψs !== nothing
                    @timeit to_d "Smoothed field ψ (Fourier)" begin
                        to_smoothed_field!(Streamfunction(), cache.longrange)
                    end
                    @timeit to_d "Interpolate to physical" begin
                        # Write interpolation output to outputs_lr[ifield + 1]
                        interpolate_to_physical!(outputs_lr[ifield += 1], cache.longrange)
                    end
                end
                if vs !== nothing
                    @timeit to_d "Smoothed field v (Fourier)" begin
                        to_smoothed_field!(Velocity(), cache.longrange)
                    end
                    @timeit to_d "Interpolate to physical" begin
                        # Write interpolation output to outputs_lr[ifield + 1]
                        interpolate_to_physical!(outputs_lr[ifield += 1], cache.longrange)
                    end
                end
                @timeit to_d "Synchronise GPU" KA.synchronize(device_lr)  # wait for the GPU to finish its work
            end
            nothing
        end
    else
        StableTasks.@spawn nothing  # empty task (returns `nothing`)
    end

    # This avoids concurrent use of `@threads :static` (two parallel sections running in
    # parallel, one for each component), which throws an error. Note that PseudoGPU is only
    # used in tests, and actually computes stuff on the CPU.
    device_lr isa PseudoGPU && wait(task_lr)

    # While the first long-range task is running, compute short-range part.
    if with_shortrange
        @timeit to "Short-range component (CPU)" begin
            @timeit to "Process point charges" process_point_charges!(cache.shortrange, pointdata)  # useful in particular for cell lists
            @timeit to "Compute Biot–Savart" add_short_range_fields!(fields, cache.shortrange, fs; LIA)
            @timeit to "Background vorticity" background_vorticity_correction!(fields, fs, params)
        end
    end

    with_longrange || return nothing

    @timeit to "Long-range component" begin
        @timeit to "Remove self-interactions (CPU)" begin
            # This is done fully on the CPU.
            if ψs !== nothing
                remove_long_range_self_interaction!(ψs, fs, Streamfunction(), params.common)
            end
            if vs !== nothing
                remove_long_range_self_interaction!(vs, fs, Velocity(), params.common)
            end
        end

        # Now wait for long-range task to finish (GPU).
        @timeit to "Wait for GPU operations" wait(task_lr)

        # Add results from long-range part.
        @timeit to "Copy output (device → host)" let
            local ifield = 0
            if ψs !== nothing
                add_long_range_output!(ψs, cache.longrange, outputs_lr[ifield += 1])
            end
            if vs !== nothing
                add_long_range_output!(vs, cache.longrange, outputs_lr[ifield += 1])
            end
        end

        foreach(KA.unsafe_free!, outputs_lr)  # manually deallocate temporary arrays
    end

    nothing
end

# Case of a list of filaments
function _compute_LIA_on_nodes!(
        fields::NamedTuple{Names, NTuple{N, V}},
        cache::BiotSavartCache,
        fs::VectorOfFilaments;
    ) where {Names, N, V <: AbstractVector{<:VectorOfVec}}
    (; to,) = cache
    (; Γ,) = cache.shortrange.params.common
    T = typeof(Γ)
    prefactor = Γ / T(4π)
    @timeit to "LIA term (only)" begin
        Threads.@threads :static for n ∈ eachindex(fs)
            f = fs[n]
            ps = @inbounds _fields_to_pairs(fields, n)
            _compute_LIA_on_nodes!(ps, cache, f; prefactor)
        end
    end
    fields
end

# Case of a single filament
function _compute_LIA_on_nodes!(
        ps::Tuple{Vararg{Pair{<:OutputField, V}}},
        cache::BiotSavartCache,
        f::AbstractFilament;
        prefactor = nothing,
    ) where {V <: VectorOfVec}
    (; params,) = cache.shortrange
    # Note: we must use the same quadrature as used when computing the globally induced terms
    (; quad, lia_segment_fraction,) = params
    (; Γ, a, Δ,) = params.common
    T = typeof(Γ)
    prefactor_ = prefactor === nothing ? (Γ / T(4π)) : prefactor
    for i ∈ eachindex(f)
        for (quantity, values) ∈ ps
            # Here `quantity` is either Velocity() or Streamfunction()
            @inbounds values[i] = local_self_induced(
                quantity, f, i, prefactor_;
                a, Δ, quad, lia_segment_fraction,
            )
        end
    end
    ps
end

end
