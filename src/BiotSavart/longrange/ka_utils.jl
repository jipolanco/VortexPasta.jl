# KernelAbstractions utils

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
