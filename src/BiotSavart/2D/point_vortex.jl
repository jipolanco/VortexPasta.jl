using LinearAlgebra: LinearAlgebra

export Vec2, PointVortices, ParamsPointVortex, SignValue, UP, DOWN, PointVortexState, ZetaHat

include("signs.jl")  # SignValue

const Vec2{T} = SVector{2, T}
const PointVortices{T} = AbstractVector{Vec2{T}}
const PointVortexVelocities{T} = AbstractVector{Vec2{T}}
const PointVortexScalars{T} = AbstractVector{T}  # usually streamfunction
const PointVortexSigns = AbstractVector{SignValue}

@doc raw"""
    ZetaHat()

Singleton type representing the vertical direction ``\hat{z}``.

In 2D, this allows to conveniently and efficiently perform operations such as
``\hat{z} \times \bm{v}`` where ``\bm{v}`` is a two-component vector (represented as a
`Vec2{T}`, i.e. an `SVector{2, T}`).

The alias `BiotSavart.ẑ = ZetaHat()` is also provided for convenience.
"""
struct ZetaHat end

const ẑ = ZetaHat()

@inline LinearAlgebra.cross(::ZetaHat, v::Vec2) = Vec2(-v[2], +v[1])
@inline LinearAlgebra.cross(v::Vec2, ::ZetaHat) = Vec2(+v[2], -v[1])

"""
    PointVortexState{T <: AbstractFloat}

Holds the instantaneous state of a point vortex system.

# Construction

    PointVortexState(points::AbstractVector{Vec2{T}}, signs::AbstractVector{SignValue})

Both vectors must have the same dimensions.

## Example

```jldoctest
julia> using VortexPasta.BiotSavart

julia> points = [Vec2(0.3, 1.54), Vec2(1.4, 5.3)]
2-element Vector{StaticArraysCore.SVector{2, Float64}}:
 [0.3, 1.54]
 [1.4, 5.3]

julia> signs = [UP, DOWN]
2-element Vector{SignValue}:
 UP::SignValue = 0x00
 DOWN::SignValue = 0x01

julia> vortices = PointVortexState(points, signs)
PointVortexState{Float64} with 2 vortices
```

# Fields

- `points::AbstractVector{Vec2{T}}`: point vortex locations;
- `signs::AbstractVector{SignValue}`: point vortex signs (`UP` or `DOWN`, see [`SignValue`](@ref)).

"""
struct PointVortexState{
        T <: AbstractFloat,
        Points <: PointVortices{T},
        Signs <: PointVortexSigns,
    }
    points::Points
    signs::Signs
    function PointVortexState(points::PointVortices{T}, signs::PointVortexSigns) where {T}
        eachindex(points) == eachindex(signs) || throw(DimensionMismatch("both vectors must have the same length"))
        new{T, typeof(points), typeof(signs)}(points, signs)
    end
end

function Base.show(io::IO, pv::PointVortexState)
    T = eltype(pv)
    N = length(pv)
    print(io, "PointVortexState{$T} with $N vortices")
end

Base.eltype(::PointVortexState{T}) where {T} = T
Base.eachindex(pv::PointVortexState) = eachindex(pv.points, pv.signs)  # fails if they have different lengths
Base.length(pv::PointVortexState) = length(eachindex(pv))
Base.isempty(pv::PointVortexState) = iszero(length(pv))
Base.copy(pv::PointVortexState) = PointVortexState(copy(pv.points), copy(pv.signs))

include("params.jl")
include("pointdata.jl")
include("pointdata_impl.jl")
include("cache.jl")
include("longrange.jl")
include("shortrange.jl")

# Define point vortex variants of 3D Biot–Savart functions.
function velocity_on_nodes(cache::PointVortexCache, vortices::PointVortexState, args...; kws...)
    vs = similar(vortices.points)
    velocity_on_nodes!(vs, cache, vortices, args...; kws...)
    vs
end

function velocity_on_nodes!(vs::PointVortexVelocities, cache::PointVortexCache, args...; kws...)
    fields = (; velocity = vs,)
    compute_on_nodes!(fields, cache, args...; kws...)
    vs
end

function background_vorticity_correction!(outputs::NamedTuple, signs::PointVortexSigns, params::ParamsCommon2D)
    # We only need to modify the streamfunction, so we do nothing if the streamfunction is
    # not present. (The velocity is not affected by the uniform background vorticity.)
    ψs = get(outputs, :streamfunction, nothing)
    ψs === nothing && return outputs
    domain_is_periodic(params) || return outputs  # correction only makes sense in periodic domain

    (; Γ, Ls, splitting,) = params

    npos = count(==(UP), signs)
    nneg = length(signs) - npos
    charge = npos - nneg
    iszero(charge) && return outputs  # nothing to do, the total charge is zero

    ω_mean = charge * Γ / prod(Ls)
    ω_back = -ω_mean  # background vorticity

    # Streamfunction associated to background vorticity
    ψ_back = ω_back * background_vorticity_correction_factor(splitting)

    AK.foreachindex(ψs) do i
        @inbounds ψs[i] = ψs[i] + ψ_back
    end

    outputs
end

function _setup_fields!(fields::NamedTuple{Names}, vortices::PointVortexState) where {Names}
    N = length(fields)
    @assert N > 0
    vs = get(fields, :velocity, nothing)
    ψs = get(fields, :streamfunction, nothing)
    vecs = filter(!isnothing, (vs, ψs))
    nfields = length(vecs)
    nfields == N || throw(ArgumentError(lazy"some of these fields were not recognised: $Names"))
    map(fields) do field
        eachindex(field) == eachindex(vortices) || throw(DimensionMismatch("wrong dimensions of output field"))
        _reset_array!(field)  # set everything to zero
        nothing
    end
end

# 2D version
function do_longrange!(
        cache::LongRangeCache, outputs::NamedTuple, pointdata_cpu::PointVortexData;
        callback_vorticity::Fvort,
    ) where {Fvort <: Function}
    (; pointdata, to) = cache  # pointdata on the device (possibly a GPU)
    TimerOutputs.reset_timer!(to)  # reset timer, since it will be merged with main timer (otherwise events will be repeated)

    # Make sure we execute this task in the GPU device chosen for long-range computations.
    # See https://cuda.juliagpu.org/dev/usage/multigpu/#Scenario-2:-Multiple-GPUs-per-process
    activate_device!(cache)

    @timeit to "Long-range component (async)" begin
        # Copy point data to the cache (possibly on a GPU).
        @assert pointdata_cpu !== pointdata  # they are different objects
        GC.@preserve pointdata begin  # see docs for KA.copyto! (it shouldn't really be needed here)
            # Resize output arrays
            noutputs = length(pointdata_cpu.points)
            map(v -> resize_no_copy!(v, noutputs), outputs)

            @timeit to "Copy point charges (host -> device)" begin
                copy_host_to_device!(pointdata.points, pointdata_cpu.points, pointdata.buf_host)
                copy_host_to_device!(pointdata.signs, pointdata_cpu.signs, pointdata.buf_host)  # we prefer transferring UPs and DOWNs (UInt8)
            end
            @assert pointdata.points !== pointdata_cpu.points  # we don't want to modify `pointdata_cpu`!
            @timeit to "Process point charges" process_point_charges!(cache)  # modifies pointdata (points and nodes)
            @timeit to "Compute charges" let
                resize_no_copy!(pointdata.charges, noutputs)
                signs_to_charges!(pointdata.charges, pointdata.signs)
            end

            # Compute vorticity in Fourier space from point data (vortex locations) -> type 1 NUFFT
            @timeit to "Vorticity to Fourier" compute_vorticity_fourier!(cache)  # reads pointdata (points and charges)

            if callback_vorticity !== identity
                @timeit to "Vorticity callback" callback_vorticity(cache)
            end

            # Interpolate streamfunction and/or velocity.
            callback_interp = get_ewald_interpolation_callback(cache)  # perform Ewald smoothing before interpolating

            # # Note: streamfunction (if enabled) must be computed before velocity.
            if hasproperty(outputs, :streamfunction)
                # Compute streamfunction from vorticity in Fourier space.
                @timeit to "Streamfunction field (Fourier)" compute_field_fourier!(Streamfunction(), cache)
                # Write interpolation output to outputs.streamfunction
                @timeit to "Interpolate to physical" interpolate_to_physical!(callback_interp, outputs.streamfunction, cache)
            end

            if hasproperty(outputs, :velocity)
                # Compute velocity from vorticity or streamfunction in Fourier space.
                @timeit to "Velocity field (Fourier)" compute_field_fourier!(Velocity(), cache)
                # Write interpolation output to outputs.velocity
                @timeit to "Interpolate to physical" interpolate_to_physical!(callback_interp, outputs.velocity, cache)
            end

            yield()  # let other tasks run (not sure if this really helps)

            # Wait for the GPU to finish its work before finishing this task.
            @timeit to "Synchronise GPU" KA.synchronize(KA.get_backend(cache))
        end
    end

    nothing
end

# 2D version
function do_shortrange!(
        cache::ShortRangeCache, outputs::NamedTuple,
        pointdata_cpu::PointVortexData,
        options::@NamedTuple{}  # this is an empty NamedTuple (for compatibility with the 3D version)
    )
    (; pointdata, to,) = cache  # pointdata on the device (possibly a GPU)
    TimerOutputs.reset_timer!(to)  # reset timer, since it will be merged with main timer (otherwise events will be repeated)

    # Make sure we execute this task in the GPU device chosen for long-range computations.
    # See https://cuda.juliagpu.org/dev/usage/multigpu/#Scenario-2:-Multiple-GPUs-per-process
    activate_device!(cache)

    @timeit to "Short-range component (async)" begin
        GC.@preserve pointdata begin  # see docs for KA.copyto! (it shouldn't really be needed here)
            # Resize output arrays
            noutputs = length(pointdata_cpu.points)
            map(v -> resize_no_copy!(v, noutputs), outputs)

            @timeit to "Copy point charges (host -> device)" begin
                copy_host_to_device!(pointdata.points, pointdata_cpu.points, pointdata.buf_host)
                copy_host_to_device!(pointdata.signs, pointdata_cpu.signs, pointdata.buf_host)  # we prefer transferring UPs and DOWNs (UInt8)
            end
            @assert pointdata.points !== pointdata_cpu.points  # we don't want to modify `pointdata_cpu`!

            @timeit to "Process point charges" process_point_charges!(cache)   # useful in particular for cell lists
            @timeit to "Fill with zeros" _reset_vectors!(outputs)  # we need to set everything to 0, since the function below _adds_ to existing values
            @timeit to "Pair interactions" add_pair_interactions_2D!(outputs, cache)
            @timeit to "Remove self-interactions" remove_self_interaction_2D!(outputs, cache)
            @timeit to "Background vorticity correction" background_vorticity_correction!(outputs, pointdata.signs, cache.params.common)

            yield()  # let other tasks run (not sure if this really helps)

            # Wait for the GPU to finish its work before finishing this task.
            @timeit to "Synchronise GPU" KA.synchronize(KA.get_backend(cache))
        end
    end

    nothing
end

function compute_on_nodes!(
        fields::NamedTuple, cache::PointVortexCache, vortices::PointVortexState;
        longrange::Bool = true, shortrange::Bool = true,
        callback_vorticity::Fvort = identity,
    ) where {Fvort <: Function}
    options = (;)  # currently point vortex has no options (passed to do_shortrange!)
    _compute_on_nodes!(
        fields, cache, vortices, options;
        with_shortrange = shortrange,
        with_longrange = longrange,
        callback_vorticity,
    )
    fields
end

# This one is called by _add_output_from_async_task! -> copy_output_values_on_nodes!
# These are the point vortex versions, which adapt to the format of the two-component velocity field
# and the one-component streamfunction field.
function _copy_output_values_on_nodes!(
        backend::GPU, op::F, vs::PointVortexVelocities{T},
        vs_d::StructVector{Vec2{T}}, buf_host::HostVector{T},
    ) where {F <: Function, T}
    for (i, us_i) in enumerate(StructArrays.components(vs_d))
        # Perform GPU -> CPU copy (device-to-host)
        @assert backend == KA.get_backend(buf_host)
        copy_device_to_host!(buf_host, us_i)  # here `us_i` is an AbstractVector{T} (usually a GPU vector)
        _copy_output_values_on_nodes!(CPU(), op, vs, i, buf_host)  # finally, perform CPU -> CPU copy
    end
    vs
end

# Copy a single component of a vector field on the CPU.
function _copy_output_values_on_nodes!(
        ::CPU, op::F, vs::PointVortexVelocities{T}, i::Int,
        vs_h::HostVector{T}, ignored = nothing,
    ) where {F <: Function, T}
    Threads.@threads for j in eachindex(vs, vs_h)
        @inbounds begin
            v⃗_prev = vs[j]
            vi_prev = v⃗_prev[i]  # old value of v⃗[i]
            vi_new = @inline op(vi_prev, vs_h[j])  # new value of v⃗[i] // typically op == +, meaning that we add to the previous value
            vs[j] = Base.setindex(v⃗_prev, vi_new, i)  # update v⃗[i]
        end
    end
end

# Case of a scalar field (usually streamfunction)
function _copy_output_values_on_nodes!(
        backend::GPU, op::F, vs::PointVortexScalars{T},
        vs_d::AbstractVector{T}, buf_host::HostVector{T},
    ) where {F <: Function, T <: Real}
    # Perform GPU -> CPU copy (device-to-host)
    @assert backend == KA.get_backend(buf_host)
    copy_device_to_host!(buf_host, vs_d)
    _copy_output_values_on_nodes!(CPU(), op, vs, buf_host)  # finally, perform CPU -> CPU copy
end

# Case of a scalar or vector field
function _copy_output_values_on_nodes!(
        ::CPU, op::F, vs::PointVortexScalars{T}, vs_h::AbstractVector{T}, ignored = nothing,
    ) where {F <: Function, T <: Union{Real, SVector}}
    Threads.@threads for j in eachindex(vs, vs_h)
        @inbounds vs[j] = @inline op(vs[j], vs_h[j])  # typically op == +, meaning that we add to the previous value
    end
end
