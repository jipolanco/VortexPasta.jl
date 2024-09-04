# KernelAbstractions utils (CPU/GPU kernels)

"""
    KernelAbstractions.get_backend(backend::LongRangeBackend) -> KernelAbstractions.Backend
    KernelAbstractions.get_backend(cache::LongRangeCache) -> KernelAbstractions.Backend

Get KernelAbstractions (KA) backend associated to a given long-range backend.

!!! note

    The word "backend" means two different things here!
    For KA, it refers to the device where kernels are executed (e.g. `CPU`, `CUDABackend`, ...).

By default this returns `KA.CPU(static = true)`, meaning that things are run on the CPU
using threads, and that a static thread assignment is used.
"""
KA.get_backend(::LongRangeBackend) = ka_default_cpu_backend()

ka_default_cpu_backend() = KA.CPU(static = true)

# Default workgroup size used for running KA kernels.
function ka_default_workgroupsize(::KA.CPU, dims::Dims)
    # On the CPU, use KA's default, which currently tries to creates blocks of 1024 work
    # items. Note that we're calling an internal function which may change in the future!
    KA.default_cpu_workgroupsize(dims)
end

function ka_default_workgroupsize(::KA.GPU, dims::Dims)
    # On the GPU, we divide the work across the first dimension and try to use 64 GPU
    # threads per workgroup.
    wgsize_wanted = 64
    x = map(one, dims)  # = (1, 1, 1) in 3D
    Base.setindex(x, min(dims[1], wgsize_wanted), 1)  # usually (wgsize_wanted, 1, 1)
end

"""
    ka_generate_kernel(kernel, backend::KA.Backend, x::AbstractArray; [workgroupsize])
    ka_generate_kernel(kernel, backend::KA.Backend, ndrange::Dims; [workgroupsize])

Generate statically sized KA kernel.

In this context, "statically sized" means that the kernel will be specifically compiled for
the dimensions of the array `x`, and will be recompiled if an array of a different size is
used later.

Here `kernel` is a KA kernel (a Julia function) annotated with the `@kernel` macro.

By default, the workgroupsize is determined automatically and may depend on the actual
backend (CPU, GPU) and on the array dimensions `ndrange = size(x)`.
"""
function ka_generate_kernel(
        kernel::F, backend::KA.Backend, ndrange::Dims;
        workgroupsize = ka_default_workgroupsize(backend, ndrange),
    ) where {F <: Function}
    kernel(backend, workgroupsize, ndrange)
end

function ka_generate_kernel(
        kernel::F, backend::KA.Backend, u::AbstractArray;
        kws...,
    ) where {F <: Function}
    ka_generate_kernel(kernel, backend, size(u); kws...)
end

## ================================================================================ ##
## This is for testing some GPU-specific code on CPUs. Used only in tests.

struct PseudoGPU <: KA.GPU end

KA.isgpu(::PseudoGPU) = false  # needed to be considered as a CPU backend by KA
KA.allocate(::PseudoGPU, args...) = KA.allocate(KA.CPU(), args...)
KA.synchronize(::PseudoGPU) = nothing
KA.copyto!(::PseudoGPU, u, v) = copyto!(u, v)
Adapt.adapt(::PseudoGPU, u::Array) = copy(u)  # simulate host → device copy (making sure arrays are not aliased)

# Convert kernel to standard CPU kernel (relies on KA internals...)
function (kernel::KA.Kernel{PseudoGPU, GroupSize, NDRange, Fun})(args...; kws...) where {GroupSize, NDRange, Fun}
    kernel_cpu = KA.Kernel{KA.CPU, GroupSize, NDRange, Fun}(KA.CPU(), kernel.f)
    kernel_cpu(args...; kws...)
end
