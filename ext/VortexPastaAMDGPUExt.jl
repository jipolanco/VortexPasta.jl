module VortexPastaAMDGPUExt

using AMDGPU: AMDGPU, HIP, ROCBackend
using VortexPasta.BiotSavart

# See:
# - VortexPastaCUDAExt.jl
# - https://github.com/JuliaGPU/AMDGPU.jl/blob/4d9bff975a46d7867b23c0d1165b05f426dd42ab/src/runtime/memory/hip.jl#L197
# - https://rocm.docs.amd.com/projects/HIP/en/docs-7.2.0/doxygen/html/group___memory.html#gab8258f051e1a1f7385f794a15300e674
function BiotSavart.pagelock!(::ROCBackend, A::DenseArray{T}) where {T}
    @debug "Pinning CPU array using AMDGPU"
    ptr = pointer(A)
    bytesize = length(A) * sizeof(T)
    flags = HIP.hipHostRegisterPortable  # "Memory is considered registered by all contexts. HIP only supports one context so this is always assumed true."
    if bytesize > 0
        HIP.hipHostRegister(ptr, bytesize, flags)
    end
    nothing
end

function BiotSavart.unpagelock!(::ROCBackend, A::DenseArray{T}) where {T}
    @debug "Unpinning CPU array using AMDGPU"
    ptr = pointer(A)
    HIP.hipHostUnregister(ptr)
    nothing
end

end
