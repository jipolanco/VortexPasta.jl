module VortexPastaOpenCLExt

using OpenCL: CLArray
using VortexPasta.BiotSavart: BiotSavart, HostVector

# unsafe_copyto! is not defined between Ptr and CLPtr, so we use unsafe_wrap instead to work
# with Array and CLArray.
function BiotSavart.copyto_ptr!(dst::CLArray, src::HostVector, n)
    @debug "Copying HostVector -> CLArray"
    GC.@preserve dst src begin
        A = unsafe_wrap(Array, pointer(src), (n,); own = false)
        copyto!(dst, A)
    end
end

function BiotSavart.copyto_ptr!(dst::HostVector, src::CLArray, n)
    @debug "Copying CLArray -> HostVector"
    GC.@preserve dst src begin
        A = unsafe_wrap(Array, pointer(dst), (n,); own = false)
        copyto!(A, src)
    end
end

end
