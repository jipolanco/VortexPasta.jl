## ========================================================================================== ##

"""
    HostVector{T, Backend <: KA.Backend} <: AbstractVector{T}

CPU vector involved in host-device transfers.

This may be used as an intermediate staging area for host-device copies.

The vector is always pagelocked (pinned) for efficient transfers.

# Construction

    HostVector{T}(undef, backend::KA.Backend, n::Int)

where:

- `T` is the element type (e.g. `Float64`);
- `backend` is a KernelAbstractions backend (such as `CUDABackend` or `ROCBackend`);
- `n` is the vector length.
"""
mutable struct HostVector{T, Backend <: KA.Backend} <: AbstractVector{T}
    data::Memory{T}  # length ≥ n
    n::Int  # logical length (exposed to the user)
    backend::Backend
end

function HostVector{T}(::UndefInitializer, backend::KA.Backend, n::Int) where {T}
    data = Memory{T}(undef, n)
    v = HostVector(data, n, backend)
    pagelock!(v)
    finalizer(unpagelock!, v)  # remove pagelock when vector is freed
    v
end

function Adapt.adapt_structure(backend::KA.Backend, v::HostVector{T}) where {T}
    HostVector{T}(undef, backend, length(v))
end

KA.get_backend(v::HostVector) = v.backend

Base.size(v::HostVector) = (v.n,)
Base.pointer(v::HostVector) = pointer(v.data)
Base.similar(v::HostVector, ::Type{S}, dims::Dims{1}) where {S} = HostVector{S}(undef, v.backend, dims[1])
Base.IndexStyle(::Type{<:HostVector}) = IndexLinear()

# These two functions should be overridden in package extensions for CUDA, AMDGPU, ...
# Note that KernelAbstractions has pagelock! but not unpagelock!, so we redefine both for consistency.
# Moreover, KA doesn't provide the flexibility we need. For example, on CUDA we want to set
# the CU_MEMHOSTREGISTER_PORTABLE flag so that memory is considered to be pinned by all
# available GPUs (useful if we're doing multi-GPU).
pagelock!(backend::KA.Backend, A::DenseArray) = nothing
unpagelock!(backend::KA.Backend, A::DenseArray) = nothing

pagelock!(v::HostVector) = pagelock!(v.backend, v.data)
unpagelock!(v::HostVector) = unpagelock!(v.backend, v.data)

# Resize HostVector, "transferring" pagelock to new memory if needed.
function resize_no_copy!(v::HostVector{T}, n) where {T}
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
        unpagelock!(v)
        v.data = data_new
        pagelock!(v)
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

function Base.copyto!(v::HostVector{T}, src::DenseArray{T}) where {T}
    n = length(src)
    resize_no_copy!(v, n)
    # This should work both when src is either a CPU or a GPU array.
    unsafe_copyto!(pointer(v), pointer(src), n)  # TODO: parallelise copy when src is on the CPU?
end

function Base.copyto!(dst::DenseArray{T}, v::HostVector{T}) where {T}
    n = length(v)
    resize_no_copy!(dst, n)
    # This should work both when src is either a CPU or a GPU array.
    unsafe_copyto!(pointer(dst), pointer(v), n)  # TODO: parallelise copy when dst is on the CPU?
end

## ========================================================================================== ##

function copy_host_to_device!(dst::AbstractVector, src::AbstractVector, buf::HostVector)
    @assert KA.get_backend(buf) == KA.get_backend(dst)
    copyto!(buf, src)  # copy host_unpinned -> host_pinned
    copyto!(dst, buf)  # copy host_pinned -> device (synchronous)
    dst
end

function copy_host_to_device!(dst::StructVector, src::StructVector, buf::HostVector)
    copy_host_to_device!(StructArrays.components(dst), StructArrays.components(src), buf)
    dst
end

function copy_host_to_device!(dst::Tuple, src::Tuple, buf::HostVector)
    foreach(dst, src) do a, b
        copy_host_to_device!(a, b, buf)
    end
    dst
end

function copy_device_to_host!(dst::AbstractVector, src::AbstractVector, buf::HostVector)
    copy_device_to_host!(buf, src)
end

function copy_device_to_host!(buf::AbstractVector, src::AbstractVector)
    @assert KA.get_backend(buf) == KA.get_backend(src)
    copyto!(buf, src)  # copy device -> host_pinned (synchronous)
    buf
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
