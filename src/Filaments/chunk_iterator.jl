public FilamentChunkIterator

# Types enabling the convenient iteration over filament nodes when working with vectors of
# filaments, especially when parallelising over multiple threads. The idea is that each
# thread deals with roughly the same number of filament _nodes_, even when the number of
# nodes per filament is very unequal across filaments.
#
# We define two iterators:
#
# - FilamentChunkIterator, which divides the list of filament nodes into chunks (1 chunk = 1 thread),
#   returning a SingleChunkIterator for each thread;
# - SingleChunkIterator, which simplifies iteration over a single chunk (by a single thread).
#
# External code only needs to know about FilamentChunkIterator.

# We iterate over filaments fs, possibly starting somewhere in the middle of filament
# fs[begin], and possibly ending somewhere in the middle of filament fs[end].
# Note that fs is a SubArray, i.e. a view of a larger vector of filaments.
struct SingleChunkIterator{Filaments <: AbstractVector{<:AbstractVector}}
    fs::Filaments         # all filaments (included those not in this chunk)
    inds::UnitRange{Int}  # indices of considered filaments
    a::Int  # index of start node in first filament (inclusive) -> iterate over nodes fs[inds[begin]][a:end]
    b::Int  # index of end node in last filament (inclusive)    -> iterate over nodes fs[inds[end]][begin:b]
end

Base.IteratorSize(::Type{<:SingleChunkIterator}) = Base.HasLength()
Base.length(it::SingleChunkIterator) = length(it.inds)  # = number of filaments included in this chunk
Base.IteratorEltype(::Type{<:SingleChunkIterator}) = Base.HasEltype()
Base.eltype(::SingleChunkIterator) = Tuple{Int, UnitRange{Int}, Int}  # = (filament_idx, node_indices, num_nodes_visited)

function Base.iterate(it::SingleChunkIterator)
    (; fs, inds, a, b) = it
    isempty(inds) && return nothing
    # First filament: start from index a
    i = first(inds)
    f = fs[i]
    if length(inds) == 1  # the chunk covers (part of) a single filament only
        node_indices = a:b
    else  # the chunk covers two or more filaments
        node_indices = a:lastindex(f)
    end::UnitRange{Int}
    # Count accumulated number of nodes until the end of the previous filament.
    num_nodes_visited = sum(j -> length(fs[j]), firstindex(fs):(i - 1); init = 0)
    # Include previously visited nodes in this filament (if a > firstindex(f)).
    num_nodes_visited += a - firstindex(f)
    ret = (i, node_indices, num_nodes_visited)::eltype(it)
    num_nodes_visited += length(node_indices)  # include nodes visited in this iteration
    state = (firstindex(inds) + 1, num_nodes_visited)  # (index of next filament, number of visited nodes)
    ret, state
end

function Base.iterate(it::SingleChunkIterator, state::Tuple)
    (; fs, inds, b) = it
    j, num_nodes_visited = state
    j == lastindex(inds) + 1 && return nothing  # we're done iterating over filaments
    i = inds[j]
    f = fs[i]
    if j == lastindex(inds)
        node_indices = firstindex(f):b
    else
        node_indices = UnitRange(eachindex(f))
    end::UnitRange{Int}
    ret = (i, node_indices, num_nodes_visited)::eltype(it)
    num_nodes_visited += length(node_indices)
    state = (j + 1, num_nodes_visited)
    ret, state
end

# ==================================================================================================== #

# Try to distribute filaments nodes over different threads so that each thread has approximately
# the same number of filament nodes (discrete points). In fact, the number of nodes per
# filament may be very unequal in practical situations, with e.g. a single filament having a
# lot of points and many other small filaments, so this can help with load balancing.

"""
    FilamentChunkIterator(
        fs::AbstractVector{<:AbstractVector};
        nchunks = Threads.nthreads(),
        full_vectors = false,
    ) -> FilamentChunkIterator

Simplifies iterating over lists of filaments, especially when parallelising over CPU threads.

This iterator divides the list of filaments into (maximum) `nchunks` chunks, so that each
chunk has roughly the same number of _nodes_ (discretisation points). The number of chunks
is by default the total number of available threads. The idea is that each thread deals with
roughly the same number of nodes, so that all of them have the same amount of work to do
(load balancing).

This iterator is particularly well adapted to the case where the number of nodes per
filament is very unequal, with e.g. a single filament having a lot of points and many other
small filaments, which is actually quite common in turbulence simulations.

To efficiently deal with that kind of case, each chunk can start or end in the middle of a
filament, meaning that a given thread may have to perform work on a subset of all nodes of a
given filament.

Optionally, if one must perform global operations on all nodes of a filament (for example
compute spline interpolation coefficients), one can pass `full_vectors = true`. In that
case, filaments will not be broken into multiple subsets, and a given filament is guaranteed
to be included in exactly a single chunk.

# Examples

Given a list of filaments `fs`, this iterator can be used as follows:

```julia
using VortexPasta.Filaments: FilamentChunkIterator

# Iterate in parallel over all filament nodes.
@sync for chunk in FilamentChunkIterator(fs)
    Threads.@spawn for (i, inds, num_nodes_visited) in chunk
        # Do some work on node indices `inds` of filament `fs[i]`.
        # `num_nodes_visited` is an integer equal to the accumulated number of nodes in
        # all previous chunks, which includes all filaments fs[begin:(i - 1)] and possibly
        # the current one, fs[i]. This may be useful in some cases.
    end
end
```

This iterator can work not only with lists of filaments, but also with any kind of vector of vectors.
For example, one could define the collection:

```julia
vs = [rand(rand(1:100)) for _ in 1:20]  # define 20 vectors of variable length (between 1 and 100 each)
@sync for chunk in FilamentChunkIterator(vs)
    # ...
end
```

Using `full_vectors = true`:

```julia
# Iterate in parallel over all filament nodes.
@sync for chunk in FilamentChunkIterator(fs; full_vectors = true)
    Threads.@spawn for (i, inds, num_nodes_visited) in chunk
        @assert inds == eachindex(fs[i])  # this is always true if full_vectors = true
        # Do some work on the whole filament `fs[i]`.
        # Similarly to before, `num_nodes_visited` is an integer equal to the accumulated
        # number of nodes in all previous chunks, which includes all filaments fs[begin:(i - 1)].
    end
end
```
"""
struct FilamentChunkIterator{Filaments <: AbstractVector{<:AbstractVector}}
    fs::Filaments
    nchunks::Int
    full_vectors::Bool
end

function FilamentChunkIterator(
        fs::AbstractVector{<:AbstractVector};
        nchunks = Threads.nthreads(),
        full_vectors::Bool = false,
    )
    FilamentChunkIterator(fs, nchunks, full_vectors)
end

Base.IteratorSize(::Type{<:FilamentChunkIterator}) = Base.SizeUnknown()  # the iterator may generate less than `nchunks` elements
Base.IteratorEltype(::Type{<:FilamentChunkIterator}) = Base.HasEltype()
Base.eltype(::FilamentChunkIterator{F}) where {F} = SingleChunkIterator{F}

function Base.iterate(it::FilamentChunkIterator)
    (; fs,) = it
    isempty(fs) && return nothing
    Np_total = sum(length, fs)  # total number of filament nodes
    Np_accumulated = zero(Np_total)
    nchunk = 0  # index of current chunk
    Base.require_one_based_indexing(fs)
    Base.require_one_based_indexing(first(fs))
    i_next = firstindex(fs)                  # first filament of next chunk (i)
    i_node_idx_next = firstindex(first(fs))  # first node of filament i_next to be considered (i_node_idx)
    state = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)
    iterate(it, state)
end

function Base.iterate(it::FilamentChunkIterator, state)
    (; fs, nchunks, full_vectors) = it
    (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated,) = state
    if nchunk == nchunks || i_next == lastindex(fs) + 1
        return nothing  # we're done iterating
    end
    checkbounds(fs[i_next], i_node_idx_next)  # i_node_idx_next is a node index of the next filament
    nchunk += 1
    # Make sure Np_accumulated_wanted is a multiple of a small power of 2. This might help
    # with false sharing issues (not sure).
    p = 16
    Np_accumulated_wanted = ((Np_total * nchunk) ÷ (nchunks * p)) * p
    if nchunk == nchunks
        Np_accumulated_wanted = Np_total  # make sure we iterate over all points by the last iteration
    end
    if Np_accumulated_wanted ≤ Np_accumulated
        # Nothing to do, iterate recursively with new value of nchunk
        state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
        return iterate(it, state_new)
    end
    i = i_next
    i_node_idx = i_node_idx_next
    j = i  # filament where this chunk ends (to be adjusted below)
    j_node_idx = i_node_idx - 1  # node index where this chunk ends (to be adjusted below)
    while Np_accumulated < Np_accumulated_wanted && j ≤ lastindex(fs)
        Np_wanted = Np_accumulated_wanted - Np_accumulated
        # Available nodes in current filament
        Np_available = lastindex(fs[j]) - j_node_idx
        if Np_available > Np_wanted
            if full_vectors
                # Stop at the end of the filament
                j_node_idx = j_node_idx + Np_available
                Np_accumulated += Np_available
                @assert Np_accumulated > Np_accumulated_wanted
                i_next = j + 1
                i_node_idx_next = firstindex(fs[j])
            else
                # Stop in the middle of the filament
                j_node_idx = j_node_idx + Np_wanted
                Np_accumulated += Np_wanted
                @assert Np_accumulated == Np_accumulated_wanted
                i_next = j  # continue on the same filament
                i_node_idx_next = j_node_idx + 1
            end
            state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
            ret = SingleChunkIterator(fs, i:j, i_node_idx, j_node_idx)
            return ret, state_new
        elseif Np_available == Np_wanted
            # Stop at the end of this filament
            j_node_idx = j_node_idx + Np_wanted
            Np_accumulated += Np_wanted
            @assert Np_accumulated == Np_accumulated_wanted
            i_next = j + 1
            i_node_idx_next = firstindex(fs[j])
            state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
            ret = SingleChunkIterator(fs, i:j, i_node_idx, j_node_idx)
            return ret, state_new
        else
            # Continue iterating over filaments
            Np_accumulated += Np_available
            j += 1  # jump to next filament
            j_node_idx = 0
        end
    end
    @assert Np_accumulated == Np_accumulated_wanted
    @assert j == lastindex(fs)
    j_node_idx = lastindex(fs[j])  # last node of last filament
    i_next = j + 1
    i_node_idx_next = firstindex(fs[j])
    state_new = (; nchunk, i_next, i_node_idx_next, Np_total, Np_accumulated)::typeof(state)
    ret = SingleChunkIterator(fs, i:j, i_node_idx, j_node_idx)
    return ret, state_new
end

# ==================================================================================================== #

function parallel_reduce(
        f::F, op::Op, fs::AbstractVector{<:AbstractVector}, nthreads;
        init = zero(number_type(fs)),
    ) where {F <: Function, Op <: Function}
    T = typeof(init)
    if nthreads == 1
        x = init
        for i in eachindex(fs)
            xnew = @inline f(i, eachindex(fs[i]))::T
            x = @inline op(x, xnew)
        end
    else
        x = _parallel_reduce(f, op, init, fs, nthreads)::T
    end
    x
end

# 1. Reduce a scalar quantity
function _parallel_reduce(f::F, op::Op, init::T, fs, nthreads) where {F, Op, T <: Threads.AtomicTypes}
    x_ref = Threads.Atomic{T}(init)
    @sync for chunk in FilamentChunkIterator(fs; nchunks = nthreads)
        Threads.@spawn let x_local = init
            for (i, inds, _) in chunk
                xnew = @inline f(i, inds)::T
                x_local = @inline op(x_local, xnew)::T
            end
            Threads.atomic_add!(x_ref, x_local)
        end
    end
    x_ref[]
end

# 2. Reduce a vector quantity (typically Vec3)
function _parallel_reduce(f::F, op::Op, init::SVector{N, T}, fs, nthreads) where {F, Op, N, T <: Threads.AtomicTypes}
    V = typeof(init)
    # Threads.Atomic only works with scalar quantities, so we create a tuple (actually SVector) of scalar Atomics.
    x_ref = map(v -> Threads.Atomic{T}(v), init)
    @sync for chunk in FilamentChunkIterator(fs; nchunks = nthreads)
        Threads.@spawn let x_local = init
            for (i, inds, _) in chunk
                xnew = @inline f(i, inds)::V
                x_local = @inline op(x_local, xnew)::V
            end
            for n in eachindex(x_ref)
                Threads.atomic_add!(x_ref[n], x_local[n])
            end
        end
    end
    map(x -> x[], x_ref)::V
end
