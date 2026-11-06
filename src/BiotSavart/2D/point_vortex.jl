export Vec2, PointVortices, ParamsPointVortex, SignValue, UP, DOWN

include("signs.jl")  # SignValue

const Vec2{T} = SVector{2, T}
const PointVortices{T} = AbstractVector{Vec2{T}}
const PointVortexVelocities{T} = AbstractVector{Vec2{T}}
const PointVortexScalars{T} = AbstractVector{T}  # usually streamfunction
const SignVector = AbstractVector{SignValue}

include("params.jl")
include("pointdata.jl")
include("pointdata_impl.jl")
include("cache.jl")
include("longrange.jl")
include("shortrange.jl")

# Define point vortex variants of 3D Biot–Savart functions.
function velocity_on_nodes(cache::PointVortexCache, vortices::PointVortices, args...; kws...)
    vs = similar(vortices)
    velocity_on_nodes!(vs, cache, vortices, args...; kws...)
    vs
end

function velocity_on_nodes!(vs::PointVortexVelocities, cache::PointVortexCache, args...; kws...)
    fields = (; velocity = vs,)
    compute_on_nodes!(fields, cache, args...; kws...)
    vs
end

function background_vorticity_correction!(fields::NamedTuple, signs::SignVector, params::ParamsPointVortex)
    # We only need to modify the streamfunction, so we do nothing if the streamfunction is
    # not present. (The velocity is not affected by the uniform background vorticity.)
    ψs = get(fields, :streamfunction, nothing)
    ψs === nothing && return fields
    domain_is_periodic(params) || return fields  # correction only makes sense in periodic domain

    (; Γ, Ls, splitting,) = params

    npos = count(==(UP), signs)
    nneg = length(signs) - npos
    charge = npos - nneg
    iszero(charge) && return fields  # nothing to do, the total charge is zero

    ω_mean = charge * Γ / prod(Ls)
    ω_back = -ω_mean  # background vorticity

    # Streamfunction associated to background vorticity
    ψ_back = ω_back * background_vorticity_correction_factor(splitting)

    @inbounds for i in eachindex(ψs)
        ψs[i] = ψs[i] + ψ_back
    end

    fields
end

function _setup_fields_2d!(fields::NamedTuple{Names}, vortices) where {Names}
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
    ) where {Fvort}
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
function do_shortrange!(cache::ShortRangeCache, outputs::NamedTuple, pointdata_cpu::PointVortexData)
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

            yield()  # let other tasks run (not sure if this really helps)

            # Wait for the GPU to finish its work before finishing this task.
            @timeit to "Synchronise GPU" KA.synchronize(KA.get_backend(cache))
        end
    end

    nothing
end

# TODO: reduce code duplication with 3D version?
function compute_on_nodes!(
        fields::NamedTuple, cache::PointVortexCache, vortices::PointVortices, signs::SignVector;
        longrange::Bool = true, shortrange::Bool = true,
        callback_vorticity::Fvort = identity,
    ) where {Fvort}
    (; to, params, pointdata) = cache
    _setup_fields_2d!(fields, vortices)

    @assert length(fields) in (1, 2)  # maximum 2 fields: streamfunction + velocity

    with_shortrange = shortrange
    with_longrange = longrange

    # This mainly copies the inputs to `pointdata`
    @timeit to "Add point charges" add_point_charges!(pointdata, vortices, signs, params)

    channel = Channel{Symbol}(2)  # 2 is the length of the channel (for :shortrange + :longrange)
    tasks = Task[]

    if with_longrange
        @debug "Computing long-range interactions asynchronously" KA.get_backend(cache.longrange)
        let cache = cache.longrange
            # Select elements of outputs with the same names as in `fields` (in this case :velocity and/or :streamfunction).
            local outputs = NamedTuple{keys(fields)}(cache.outputs)
            # Compute long-range part asynchronously (e.g. on a GPU).
            local task = Threads.@spawn :interactive try
                do_longrange!(cache, outputs, pointdata; callback_vorticity)
            finally
                put!(channel, :longrange)  # this notifies that the task has completed (or failed)
            end
            KA.get_backend(cache) isa CPU && wait(task)  # this helps get more accurate timings
            push!(tasks, task)
        end
    end

    if with_shortrange
        @debug "Computing short-range interactions asynchronously" KA.get_backend(cache.shortrange)
        let cache = cache.shortrange
            # Select elements of outputs with the same names as in `fields` (in this case :velocity and/or :streamfunction).
            local outputs = NamedTuple{keys(fields)}(cache.outputs)
            # Compute short-range part asynchronously (e.g. on a GPU).
            local task = Threads.@spawn :interactive try
                do_shortrange!(cache, outputs, pointdata)
            finally
                put!(channel, :shortrange)  # this notifies that the task has completed (or failed)
            end
            KA.get_backend(cache) isa CPU && wait(task)  # this helps get more accurate timings
            push!(tasks, task)
        end
    end

    # Perform CPU-only operations associated to the short-range part (this differentiation is kind of arbitrary).
    if with_shortrange
        @debug "Computing synchronous operations"
        @timeit to "CPU-only operations (synchronous)" begin
            # Perform other CPU-only operations which we include in the short-range part (this choice is kind of arbitrary).
            background_vorticity_correction!(fields, signs, params)
        end
    end

    # Now wait for asynchronous long-range and short-range operations to finish.
    tree_point_base = [t.name for t in to.timer_stack]  # to merge timers (https://github.com/KristofferC/TimerOutputs.jl/issues/143)
    ntasks = length(tasks)
    while ntasks > 0
        @timeit to "Wait for async task to finish" begin
            taskname = take!(channel)::Symbol  # wait for first async task to finish
            @debug "Task finished: $taskname"
            ntasks -= 1
        end
        # Add results from asynchronous task.
        @timeit to "Copy output (device -> host)" let
            if taskname == :shortrange
                timer = cache.shortrange.to
                _add_output_from_async_task!(fields, cache.shortrange.outputs, cache.shortrange.pointdata.buf_host)
            elseif taskname == :longrange
                timer = cache.longrange.to
                _add_output_from_async_task!(fields, cache.longrange.outputs, cache.longrange.pointdata.buf_host)
            end
            # Add timings from asynchronous computations.
            TimerOutputs.merge!(to, timer; tree_point = tree_point_base)  # https://github.com/KristofferC/TimerOutputs.jl/issues/143
        end
    end

    foreach(fetch, tasks)  # if a task failed with an error, this will rethrow the error

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
    ) where {F <: Function, T}
    # Perform GPU -> CPU copy (device-to-host)
    @assert backend == KA.get_backend(buf_host)
    copy_device_to_host!(buf_host, vs_d)
    _copy_output_values_on_nodes!(CPU(), op, vs, buf_host)  # finally, perform CPU -> CPU copy
end

# Case of a scalar field (usually streamfunction)
function _copy_output_values_on_nodes!(
        ::CPU, op::F, vs::PointVortexScalars{T}, vs_h::AbstractVector{T}, ignored = nothing,
    ) where {F <: Function, T}
    Threads.@threads for j in eachindex(vs, vs_h)
        @inbounds vs[j] = @inline op(vs[j], vs_h[j])  # typically op == +, meaning that we add to the previous value
    end
end
