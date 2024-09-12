module VortexPastaThreadPinningExt

using VortexPasta.BiotSavart: BiotSavart, FINUFFTBackend
using ThreadPinning: ThreadPinning

# Returns true if threads are currently pinned.
threads_are_pinned(tids) = all(i -> ThreadPinning.ispinned(threadid = i), tids)

# Note: we only redefine this function for FINUFFTBackend (but not for CuFINUFFTBackend,
# which doesn't use CPU threads).
# This is a workaround for apparent conflicts between OpenMP (used by FINUFFT) and
# ThreadPinning.pinthreads.
# See also the implementation of ThreadPinning.with_pinthreads (which basically does the opposite).
function BiotSavart.finufft_unpin_threads(f::F, ::FINUFFTBackend) where {F}
    threadpool = :default
    tids = ThreadPinning.threadids(; threadpool)
    was_pinned = threads_are_pinned(tids)   # threads are currently pinned?
    if was_pinned
        # This is based on the ThreadPinning.with_pinthreads implementation.
        masks_prior = [ThreadPinning.getaffinity(; threadid = i) for i ∈ tids]  # for restoring affinity mask later (not sure if needed, but just in case)
        cpuids_prior = ThreadPinning.getcpuids()  # for restoring pinning later
        ThreadPinning.unpinthreads(; threadpool)  # disable pinning before calling f
    end
    try
        f()
    finally
        if was_pinned
            ThreadPinning.pinthreads(cpuids_prior)  # restore pinning
            for (i, threadid) ∈ pairs(tids)
                ThreadPinning.setaffinity(masks_prior[i]; threadid)  # restore affinity mask (not sure if needed)
            end
            @assert threads_are_pinned(tids)
        end
    end
end

end
