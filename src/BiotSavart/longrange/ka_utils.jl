# KernelAbstractions utils (CPU/GPU kernels)

"""
    KernelAbstractions.get_backend(backend::LongRangeBackend) -> KernelAbstractions.Backend
    KernelAbstractions.get_backend(cache::LongRangeCache) -> KernelAbstractions.Backend

Get KernelAbstractions (KA) backend associated to a given long-range backend.

The word "backend" means two different things here! For KA, it refers to the device where
kernels are executed (e.g. `CPU`, `CUDABackend`, ...).

By default this returns `KA.CPU(static = true)`, meaning that things are run on the CPU
using threads, and that a static thread assignment is used.
"""
KA.get_backend(::LongRangeBackend) = KA.CPU(static = true)

# Default workgroup size used for running KA kernels.
# This may be changed in the future (or overridden in a package extension) so that certain
# backends (e.g. GPUs) use a different workgroup size.
function ka_default_workgroupsize(::KA.Backend, dims::Dims)
    wgsize_wanted = 64  # TODO: adapt for CPU/GPU?
    groupsize = map(one, dims)  # = (1, 1, 1) in 3D
    Base.setindex(groupsize, min(groupsize[1], wgsize_wanted), 1)  # usually (wgsize_wanted, 1, 1)
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
