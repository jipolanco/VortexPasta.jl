## ========================================================================================== ##

"""
    HostVector{T} <: AbstractVector{T}

CPU vector involved in host-device transfers.

This may be used as an intermediate staging area for host-device copies.

The vector is always pagelocked (pinned) for efficient transfers.
"""
mutable struct HostVector{T} <: AbstractVector{T}
    data::Memory{T}  # length ≥ n
    n::Int  # logical length (exposed to the user)
end

function allocate(backend::KA.Backend, ::Type{HostVector{T}}, n::Int) where {T}
    data = Memory{T}(undef, n)
    v = HostVector(data, n)
    pagelock!(backend, v)
    finalizer(x -> unpagelock!(backend, x), v)  # remove pagelock when vector is freed
    v
end

Base.size(v::HostVector) = (v.n,)
Base.pointer(v::HostVector) = pointer(v.data)
Base.similar(v::HostVector, ::Type{S}, dims::Dims{1}) where {S} = HostVector(similar(v.data, S, dims), dims[1], false)
Base.IndexStyle(::Type{<:HostVector}) = IndexLinear()

# These two functions should be overridden in package extensions for CUDA, AMDGPU, ...
# Note that KernelAbstractions has pagelock! but not unpagelock!, so we redefine both for consistency.
# Moreover, KA doesn't provide the flexibility we need. For example, on CUDA we want to set
# the CU_MEMHOSTREGISTER_PORTABLE flag so that memory is considered to be pinned by all
# available GPUs (useful if we're doing multi-GPU).
pagelock!(backend::KA.Backend, A::DenseArray) = nothing
unpagelock!(backend::KA.Backend, A::DenseArray) = nothing

pagelock!(backend::KA.Backend, v::HostVector) = pagelock!(backend, v.data)
unpagelock!(backend::KA.Backend, v::HostVector) = unpagelock!(backend, v.data)

# Resize HostVector, "transferring" pagelock to new memory if needed.
function resize_no_copy!(backend::KA.Backend, v::HostVector{T}, n) where {T}
    n == length(v) && return v  # nothing to do
    capacity = length(v.data)
    @assert v.n ≤ capacity  # this should always be the case
    if n > capacity
        # Reallocate data (dropping old data).
        # Note that old data is not copied (we assume we don't need it).
        capacity_new = nextpow(2, n)
        @assert capacity_new ≥ n
        data_new = Memory{T}(undef, capacity_new)
        # We need to redo the pagelock since the pointer and data size have changed.
        unpagelock!(backend, v)
        v.data = data_new
        pagelock!(backend, v)
    end
    v.n = n  # change the logical vector size
    @assert n ≤ length(v.data)
    v
end

Base.@propagate_inbounds function Base.getindex(v::HostVector, i::Int)
    @boundscheck checkbounds(v, i)
    v.data[i]
end

Base.@propagate_inbounds function Base.setindex!(v::HostVector, x, i::Int)
    @boundscheck checkbounds(v, i)
    v.data[i] = x
end

## ========================================================================================== ##

function copy_host_to_device!(dst::AbstractVector, src::HostVector)
    resize_no_copy!(dst, length(src))
    # backend = KA.get_backend(dst)
    # KA.copyto!(backend, dst, src)  # async
    copyto!(dst, src)  # sync
    dst
end

function copy_host_to_device!(dst::StructVector, src::StructVector)
    copy_host_to_device!(StructArrays.components(dst), StructArrays.components(src))
    dst
end

function copy_host_to_device!(dst::Tuple, src::Tuple)
    foreach(dst, src) do a, b
        copy_host_to_device!(a, b)
    end
    dst
end

function copy_device_to_host!(dst::HostVector, src::AbstractVector)
    backend = KA.get_backend(src)
    resize_no_copy!(backend, dst, length(src))
    # KA.copyto!(backend, dst, src)  # asynchronous copy (on CUDA at least)
    copyto!(dst, src)  # synchronous copy
    dst
end

function copy_device_to_host!(dst::StructVector, src::StructVector)
    copy_device_to_host!(StructArrays.components(dst), StructArrays.components(src))
    dst
end

function copy_device_to_host!(dst::Tuple, src::Tuple)
    foreach(dst, src) do a, b
        copy_device_to_host!(a, b)
    end
    dst
end
