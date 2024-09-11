module VortexPastaThreadPinningExt

using VortexPasta.BiotSavart: BiotSavart, FINUFFTBackend
using ThreadPinning: ThreadPinning

# Note: we only redefine this function for FINUFFTBackend (but not for CuFINUFFTBackend,
# which doesn't use CPU threads).
# This is a workaround for apparent conflicts between OpenMP (used by FINUFFT) and
# ThreadPinning.pinthreads.
function BiotSavart.finufft_unpin_threads(f::F, ::FINUFFTBackend) where {F}
    was_pinned = ThreadPinning.ispinned()   # threads are currently pinned?
    if was_pinned
        cpuids = ThreadPinning.getcpuids()  # for restoring pinning later
        ThreadPinning.unpinthreads()        # disable pinning before calling f
    end
    try
        f()
    finally
        if was_pinned
            ThreadPinning.pinthreads(cpuids)  # restore pinning
            @assert ThreadPinning.ispinned()
        end
    end
end

end
