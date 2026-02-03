module VortexPastaCUDAExt

using CUDA: CUDA, CUDABackend
using VortexPasta.BiotSavart

# Adapted from https://github.com/JuliaGPU/CUDA.jl/blob/e2eab84e89085b5298d44634858993939479f0ac/lib/cudadrv/memory.jl#L167
# See also https://docs.nvidia.com/cuda/cuda-driver-api/group__CUDA__MEM.html#group__CUDA__MEM_1gf0a9fe11544326dabd743b7aa6b54223
function BiotSavart.pagelock!(::CUDABackend, A::DenseArray{T}) where {T}
    @debug "Pinning CPU array using CUDA"
    ptr = pointer(A)
    bytesize = length(A) * sizeof(T)
    flags = CUDA.MEMHOSTREGISTER_PORTABLE  # "The memory returned by this call will be considered as pinned memory by all CUDA contexts, not just the one that performed the allocation."
    if !isempty(A)
        CUDA.cuMemHostRegister_v2(ptr, bytesize, flags)
    end
    nothing
end

function BiotSavart.unpagelock!(::CUDABackend, A::DenseArray{T}) where {T}
    @debug "Unpinning CPU array using CUDA"
    ptr = pointer(A)
    if !isempty(A)
        CUDA.cuMemHostUnregister(ptr)
    end
    nothing
end

end
