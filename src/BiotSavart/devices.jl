# Tools for working with CPU and GPU devices (via KernelAbstractions.jl)

"""
    KernelAbstractions.get_backend(backend::AbstractBackend) -> Union{CPU, KA.Backend}
    KernelAbstractions.get_backend(cache::ShortRangeCache) -> Union{CPU, KA.Backend}
    KernelAbstractions.get_backend(cache::LongRangeCache) -> Union{CPU, KA.Backend}

Get device type (CPU or KernelAbstractions GPU backend) associated to a given short-range or long-range backend.

!!! note

    The word "backend" means two different things here!
    For KA, it refers to the type of device where kernels are executed (e.g. `CUDABackend`, `POCLBackend`, ...).

By default this returns `BiotSavart.CPU()`, meaning that things are run on the CPU using threads.
"""
KA.get_backend(::AbstractBackend) = CPU()

"""
    KernelAbstractions.device(backend::AbstractBackend) -> Int
    KernelAbstractions.device(cache::ShortRangeCache) -> Int
    KernelAbstractions.device(cache::LongRangeCache) -> Int

Return the device id (in `1:ndevices`) where short-range or long-range computations are run.

This can make sense when running on GPUs, where one may want to take advantage of multiple
available GPUs on the same machine. On CPUs the device id is generally `1`.
"""
KA.device(::AbstractBackend) = 1

"""
    BiotSavart.activate_device!(backend::AbstractBackend)
    BiotSavart.activate_device!(cache::ShortRangeCache)
    BiotSavart.activate_device!(cache::LongRangeCache)

Activate KernelAbstractions device associated to backend or cache.
"""
function activate_device!(cache)
    backend = KA.get_backend(cache)
    device = KA.device(cache)
    KA.device!(backend, device)
end

@inline function default_workgroupsize(ndrange::Dims, max_items::Int = 512)
    n = min(max_items, ndrange[1])
    (n, default_workgroupsize(Base.tail(ndrange), max_items ÷ n)...)
end

@inline default_workgroupsize(::Tuple{}, max_items::Int) = ()

# Default workgroup size used for running KA kernels.
function ka_default_workgroupsize(::DeviceType, ndrange::Dims)
    # Try to use up to 512 work items (GPU threads) per workgroup.
    default_workgroupsize(ndrange, 512)
end

"""
    ka_generate_kernel(kernel, backend::KA.Backend, x::AbstractArray; [workgroupsize])

Generate statically sized KA kernel.

In this context, "statically sized" means that the kernel will be specifically compiled for
the dimensions of the array `x`, and will be recompiled if an array of a different size is
used later.

Here `kernel` is a KA kernel (a Julia function) annotated with the `@kernel` macro.

By default, the workgroupsize is determined automatically and may depend on the actual
backend (CPU, GPU) and on the array dimensions `ndrange = size(x)`.
"""
function ka_generate_kernel(
        kernel::F, backend::DeviceType, u::AbstractArray;
        workgroupsize = ka_default_workgroupsize(backend, size(u)),
    ) where {F <: Function}
    ndrange = size(u)
    backend_ka = Devices.to_ka_backend(backend)  # converts CPU() -> POCLBackend()
    kernel(backend_ka, workgroupsize, ndrange)
end

# Resize vector trying to avoid copy when N is larger than the original length.
# In other words, we completely discard the original contents of x, which is not the
# original behaviour of resize!. This can save us some device-to-device copies.
# Function copied from NonuniformFFTs.jl.
function resize_no_copy!(x, N)
    if length(x) ≠ N
        resize!(x, 0)
        resize!(x, N)
    end
    x
end
