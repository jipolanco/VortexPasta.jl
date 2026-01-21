# KernelAbstractions utils (CPU/GPU kernels)

"""
    KernelAbstractions.get_backend(backend::AbstractBackend) -> KernelAbstractions.Backend
    KernelAbstractions.get_backend(cache::ShortRangeCache) -> KernelAbstractions.Backend
    KernelAbstractions.get_backend(cache::LongRangeCache) -> KernelAbstractions.Backend

Get KernelAbstractions (KA) backend associated to a given short-range or long-range backend.

!!! note

    The word "backend" means two different things here!
    For KA, it refers to the type of device where kernels are executed (e.g. `CPU`, `CUDABackend`, ...).

By default this returns `KA.CPU()`, meaning that things are run on the CPU using threads.
"""
KA.get_backend(::AbstractBackend) = ka_default_cpu_backend()

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

ka_default_cpu_backend() = KA.CPU()

# Default workgroup size used for running KA kernels.
function ka_default_workgroupsize(::KA.CPU, dims::Dims)
    # On the CPU, use KA's default, which currently tries to creates blocks of 1024 work
    # items. Note that we're calling an internal function which may change in the future!
    KA.default_cpu_workgroupsize(dims)
end

function ka_default_workgroupsize(::KA.GPU, dims::Dims)
    # On the GPU, we divide the work across the first dimension and try to use 512 GPU
    # threads per workgroup.
    wgsize_wanted = 512
    x = map(one, dims)  # = (1, 1, 1) in 3D
    Base.setindex(x, min(dims[1], wgsize_wanted), 1)  # usually (wgsize_wanted, 1, 1)
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
        kernel::F, backend::KA.Backend, u::AbstractArray;
        workgroupsize = ka_default_workgroupsize(backend, size(u)),
    ) where {F <: Function}
    ndrange = size(u)
    kernel(backend, workgroupsize, ndrange)
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

## ================================================================================ ##
## This is for testing some GPU-specific code on CPUs. Used only in tests.

struct PseudoGPU <: KA.GPU end

KA.isgpu(::PseudoGPU) = false  # needed to be considered as a CPU backend by KA
KA.allocate(::PseudoGPU, ::Type{T}, dims::Dims) where {T} = KA.allocate(KA.CPU(), T, dims)
KA.synchronize(::PseudoGPU) = nothing

# Emulate multiple devices using task-local storage
KA.ndevices(::PseudoGPU) = 4
KA.device(::PseudoGPU) = get(task_local_storage(), :PseudoGPU_device, 1)::Int
function KA.device!(backend::PseudoGPU, device::Int)
    1 ≤ device ≤ KA.ndevices(backend) || throw(ArgumentError("invalid device"))
    task_local_storage(:PseudoGPU_device, device)
end

KA.copyto!(::PseudoGPU, u, v) = copyto!(u, v)
Adapt.adapt(::PseudoGPU, u::Array) = copy(u)  # simulate host → device copy (making sure arrays are not aliased)

# Convert kernel to standard CPU kernel (relies on KA internals...)
function (kernel::KA.Kernel{PseudoGPU, GroupSize, NDRange, Fun})(args...; kws...) where {GroupSize, NDRange, Fun}
    kernel_cpu = KA.Kernel{KA.CPU, GroupSize, NDRange, Fun}(KA.CPU(), kernel.f)
    kernel_cpu(args...; kws...)
end
