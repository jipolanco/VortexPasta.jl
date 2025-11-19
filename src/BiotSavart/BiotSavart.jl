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
    velocity_on_nodes, velocity_on_nodes!,
    compute_on_nodes!,
    CPU,  # from KernelAbstractions
    reset_timer!  # from TimerOutputs

using ..Constants: Zero, Infinity, ∞

using ..Quadratures:
    Quadratures, quadrature, NoQuadrature, GaussLegendre, AdaptiveTanhSinh,
    AbstractQuadrature, StaticSizeQuadrature, PreallocatedQuadrature

using ..Filaments:
    Filaments, AbstractFilament, ClosedFilament, Segment, CurvatureBinormal,
    Vec3, Derivative,
    knots, nodes, segments, integrate

using Adapt: Adapt, adapt

using KernelAbstractions:
    KernelAbstractions,  # importing this avoids docs failure
    KernelAbstractions as KA, @kernel, @index, @Const,
    CPU, GPU

using StableTasks: StableTasks

using Bumper: Bumper, @no_escape, @alloc
using StructArrays: StructArrays, StructVector, StructArray
using TimerOutputs: TimerOutputs, TimerOutput, @timeit, reset_timer!

abstract type OutputField end
struct Streamfunction <: OutputField end
struct Velocity <: OutputField end

"""
    AbstractBackend

Denotes a "backend" for short-range or long-range Ewald computations.

See [`ShortRangeBackend`](@ref) and [`LongRangeBackend`](@ref) for more details.
"""
abstract type AbstractBackend end

const VectorOfFilaments = AbstractVector{<:AbstractFilament}
const VectorOfVec = AbstractVector{<:Vec3}
const VectorOfPositions = VectorOfVec
const VectorOfVelocities = VectorOfVec
const AllFilamentVelocities = AbstractVector{<:VectorOfVelocities}

include("ka_utils.jl")
include("pointdata.jl")
include("types_shortrange.jl")
include("types_longrange.jl")
include("params.jl")
include("cache.jl")
include("autotune.jl")

include("shortrange/shortrange.jl")
include("longrange/longrange.jl")

"""
    velocity_on_nodes(cache::BiotSavartCache, fs::AbstractVector{<:AbstractFilament}; kws...) -> vs

Compute velocity induced by vortex filaments on filament nodes.

Returns a `vs` vector containing the velocities on filament nodes.

See [`velocity_on_nodes!`](@ref) for more details.
"""
function velocity_on_nodes(
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament};
        kws...,
    )
    vs = map(similar ∘ nodes, fs)
    velocity_on_nodes!(vs, cache, fs; kws...)
    vs
end

"""
    velocity_on_nodes!(
        vs::AbstractVector{<:AbstractVector{<:Vec3}},
        cache::BiotSavartCache,
        fs::AbstractVector{<:AbstractFilament};
        kws...,
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

See [`compute_on_nodes!`](@ref) for a list of accepted keyword arguments.
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
    nothing
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
⟨ \bm{ω} ⟩ = \frac{Γ}{V} ∮_{\mathcal{C}} \mathrm{d}\bm{s},
```

where $V = L^3$ is the periodic domain volume.

This quantity is always zero when all filaments are closed. It can be non-zero if there are
infinite unclosed filaments which are not compensated by oppositely-oriented infinite filaments.

When the mean filament vorticity is non-zero, one must compensate it by adding a uniform background
vorticity ``\bm{ω}_{0} = -⟨ \bm{ω} ⟩``, so that the total circulation around the
volume is zero.
This should be interpreted by considering that computations are performed in a **rotating
reference frame** with rotation rate ``\bm{Ω} = 2 ⟨ \bm{ω} ⟩ = -2 \bm{ω}_{0}``.

In practice, in the long-range component of Ewald summation, this compensation is done
implicitly by setting the zero mode ``\hat{\bm{\omega}}(\bm{k} = \bm{0}) = \bm{0}``. As for
the short-range component, only the streamfunction needs to be corrected, while the velocity
is not affected by the background vorticity.

The correction to the short-range streamfunction is given by

```math
\bm{ψ}^{<}_{0}
= \frac{\bm{ω}_{0}}{4π} ∫_{\mathbb{R}^3} \frac{\operatorname{erfc}(α|\bm{r}|)}{|\bm{r}|} \, \mathrm{d}^3\bm{r}
= \frac{\bm{ω}_{0}}{4 α^2}
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
    p⃗_total = sum(Filaments.end_to_end_offset, fs)  # p⃗_total[i] is expected to be a multiple of Ls[i]
    if p⃗_total + Vec3(Ls) ≈ Vec3(Ls)  # this means p⃗_total is basically 0
        return fields  # there's nothing to do if the mean vorticity is zero
    end

    ω⃗_mean = p⃗_total * (Γ / prod(Ls))
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
Moreover, one can call [`get_longrange_field_fourier`](@ref) to get the vorticity field from
within the callback. Of course, this callback is never be called if long-range computations are disabled.

When the callback is called, the available field is simply a Fourier-truncated vorticity, as
the coarse-graining due to Ewald's method has not been applied yet.

An example of how to compute the (large-scale) kinetic energy associated to the
Fourier-truncated vorticity field:

```julia
using Adapt: adapt  # useful in case FFTs are computed on the GPU

E_from_vorticity = Ref(0.0)  # "global" variable updated when calling compute_on_nodes!

function callback_vorticity(cache::LongRangeCache)
    (; field, wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache)
    @assert state.quantity == :vorticity
    @assert state.smoothing_scale == 0  # unsmoothed field
    # To make things simple, we copy data to the CPU if it's on the GPU.
    wavenumbers = adapt(Array, wavenumbers)
    uhat = adapt(Array, field)::NTuple{3}  # (ωx, ωy, ωz) in Fourier space
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)  # this depends on the long-range backend
    @assert with_hermitian_symmetry == wavenumbers[1][end] > 0
    E = 0.0
    for I ∈ CartesianIndices(uhat[1])
        k⃗ = map(getindex, wavenumbers, Tuple(I))
        kx = k⃗[1]
        factor = (!with_hermitian_symmetry || kx == 0) ? 0.5 : 1.0
        k² = sum(abs2, k⃗)
        if !iszero(k²)
            ω⃗ = (uhat[1][I], uhat[2][I], uhat[3][I])  # Fourier coefficient of the vorticity
            E += factor * sum(abs2, ω⃗) / k²
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
    _compute_on_nodes!(fields, cache, fs; LIA, kws...)
    fields
end

function do_longrange!(
        cache::LongRangeCache, outputs::NamedTuple, pointdata_cpu;
        callback_vorticity::Fvort,
    ) where {Fvort}
    (; pointdata_d, to,) = cache  # pointdata on the device (possibly a GPU)

    TimerOutputs.reset_timer!(to)  # reset timer, since it will be merged with main timer (otherwise events will be repeated)

    # Make sure we execute this task in the GPU device chosen for long-range computations.
    # See https://cuda.juliagpu.org/dev/usage/multigpu/#Scenario-2:-Multiple-GPUs-per-process
    ka_backend = KA.get_backend(cache)  # KA backend used for long-range computations (e.g. CUDABackend)
    device_id = KA.device(cache)        # in 1:ndevices
    KA.device!(ka_backend, device_id)   # set the device

    @timeit to "Long-range component (async)" begin
        # Copy point data to the cache (possibly on a GPU).
        @assert pointdata_cpu !== pointdata_d  # they are different objects
        @timeit to "Copy point charges (host → device)" begin
            copy!(pointdata_d, pointdata_cpu)  # H2D copy
        end
        @timeit to "Process point charges" begin
            process_point_charges!(cache)  # modifies pointdata_d (points and nodes)
        end

        # Compute vorticity in Fourier space from point data (vortex locations) -> type 1 NUFFT
        @timeit to "Vorticity to Fourier" begin
            compute_vorticity_fourier!(cache)  # reads pointdata_d (points and charges)
        end
        if callback_vorticity !== identity
            @timeit to "Vorticity callback" begin
                callback_vorticity(cache)
            end
        end

        # Interpolate streamfunction and/or velocity.
        callback_interp = get_ewald_interpolation_callback(cache)  # perform Ewald smoothing before interpolating

        if haskey(outputs, :streamfunction)
            @timeit to "Streamfunction field (Fourier)" begin
                # Compute streamfunction from vorticity in Fourier space.
                compute_field_fourier!(Streamfunction(), cache)
            end
            @timeit to "Interpolate to physical" begin
                # Write interpolation output to outputs.streamfunction
                interpolate_to_physical!(callback_interp, outputs.streamfunction, cache)
            end
        end

        if haskey(outputs, :velocity)
            @timeit to "Velocity field (Fourier)" begin
                # Compute velocity from vorticity or streamfunction in Fourier space.
                compute_field_fourier!(Velocity(), cache)
            end
            @timeit to "Interpolate to physical" begin
                # Write interpolation output to outputs.velocity
                interpolate_to_physical!(callback_interp, outputs.velocity, cache)
            end
        end

        @timeit to "Synchronise GPU" KA.synchronize(ka_backend)  # wait for the GPU to finish its work
    end

    nothing
end

# CPU/GPU version
# We compute short-range (CPU) and long-range (GPU) asynchronously, so that both components
# work at the same time.
function _compute_on_nodes!(
        fields::NamedTuple, cache, fs;
        LIA = Val(true),
        longrange = true,
        shortrange = true,
        callback_vorticity::Fvort = identity,
    ) where {Fvort}
    (; to, params, pointdata,) = cache
    (; quad, Ls,) = params
    _setup_fields!(fields, fs)

    @assert length(fields) ∈ (1, 2)  # maximum 2 fields: streamfunction + velocity

    with_shortrange = shortrange
    with_longrange = longrange && cache.longrange !== NullLongRangeCache()

    # This is used by short and long-range computations.
    @timeit to "Add point charges" add_point_charges!(pointdata, fs, Ls, quad)  # done on the CPU

    # Allocate temporary arrays on the GPU for interpolation outputs (manually deallocated later).
    if with_longrange
        noutputs = sum(length, fs)  # total number of interpolation points
        # Select elements of outputs_d with the same names as in `fields` (in this case :velocity and/or :streamfunction).
        outputs_lr = NamedTuple{keys(fields)}(cache.longrange.outputs_d)
        foreach(v -> resize_no_copy!(v, noutputs), outputs_lr)
    end

    # Compute long-range part asynchronously on the GPU.
    task_lr = if with_longrange
        StableTasks.@spawn begin
            do_longrange!(cache.longrange, outputs_lr, pointdata; callback_vorticity)::Nothing  # should return nothing for type stability
        end
    else
        StableTasks.@spawn nothing  # empty task (returns `nothing`)
    end

    # While the first long-range task is running, compute short-range part.
    if with_shortrange
        @timeit to "Short-range component" begin
            @timeit to "Process point charges" process_point_charges!(cache.shortrange, pointdata)  # useful in particular for cell lists
            @timeit to "Compute Biot–Savart" add_short_range_fields!(fields, cache.shortrange, fs; LIA)
            @timeit to "Background vorticity" background_vorticity_correction!(fields, fs, params)
            @timeit to "Remove self-interactions (CPU)" begin
                # This is done fully on the CPU.
                if hasproperty(fields, :streamfunction)
                    remove_self_interaction!(fields.streamfunction, fs, Streamfunction(), params.common)
                end
                if hasproperty(fields, :velocity)
                    remove_self_interaction!(fields.velocity, fs, Velocity(), params.common)
                end
            end

        end
    end

    with_longrange || return nothing

    @timeit to "Long-range component" begin
        # Wait for long-range task to finish (GPU).
        @timeit to "Wait for async operations" wait(task_lr)

        # Add results from long-range part.
        @timeit to "Copy output (device → host)" let
            if hasproperty(fields, :streamfunction)
                copy_long_range_output!(+, fields.streamfunction, cache.longrange, outputs_lr.streamfunction)
            end
            if hasproperty(fields, :velocity)
                copy_long_range_output!(+, fields.velocity, cache.longrange, outputs_lr.velocity)
            end
        end

        TimerOutputs.merge!(to, cache.longrange.to; tree_point = ["Long-range component"])
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
        Threads.@threads :dynamic for n ∈ eachindex(fs)
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
    Threads.@threads :dynamic for i ∈ eachindex(f)
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
