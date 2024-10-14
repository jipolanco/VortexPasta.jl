module VortexPastaAMDGPUExt

using VortexPasta.BiotSavart
using AMDGPU: ROCBackend, ROCArray
using KernelAbstractions: KernelAbstractions as KA

## ========================================================================================== ##

# Workaround KA.copyto! currently not working with host arrays allocated via Bumper.jl.

function BiotSavart.copyto_bumper!(backend::ROCBackend, dst::AbstractArray, src::AbstractArray)
    # KA.copyto!(backend, dst, src_bumper)  # works on CPU and CUDA, but currently not on AMDGPU
    _copyto_bumper!(dst, src)
end

# Host-to-device copy
function _copyto_bumper!(dst::ROCArray, src::AbstractArray)
    src_pin = unsafe_wrap(ROCArray, pointer(src), size(src))  # this should pin memory, enabling asynchronous copies
    GC.@preserve src src_pin begin
        # Since AMDGPU.jl currently defines copyto! for standard Array types, we "convert"
        # the input CPU array (usually an UnsafeArray) onto a standard Array.
        u = unsafe_wrap(Array, pointer(src), size(src))
        copyto!(dst, 1, u, 1, length(dst))
    end
    dst
end

# Device-to-host copy
function _copyto_bumper!(dst::AbstractArray, src::ROCArray)
    # We use async = true to be consistent with the CUDA implementation of KA.copyto!.
    # Currently, in AMDGPU.jl, the `async` keyword argument is only accepted for device-to-host
    # copies (in other cases copies are *always* asynchronous).
    dst_pin = unsafe_wrap(ROCArray, pointer(dst), size(dst))  # this should pin memory, enabling asynchronous copies
    GC.@preserve dst dst_pin begin
        # Since AMDGPU.jl currently defines copyto! for standard Array types, we "convert"
        # the input CPU array (usually an UnsafeArray) onto a standard Array.
        v = unsafe_wrap(Array, pointer(dst), size(dst))
        copyto!(v, 1, src, 1, length(dst); async = true)
    end
    dst
end

## ========================================================================================== ##

end
