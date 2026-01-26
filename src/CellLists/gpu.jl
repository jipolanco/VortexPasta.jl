using KernelAbstractions: @kernel, @index, @Const

@inline function determine_cell_index_gpu(x, rcut, L, N)
    # (1) Make sure x is in [0, L).
    # This is based on NonuniformFFTs.to_unit_cell_gpu.
    r = rem(x, L)   # note: rem(x, y) translates to fmodf/fmod on CUDA (see https://github.com/JuliaGPU/CUDA.jl/blob/4067511b2b472be9fb30164d1ed23caa354c1fcb/src/device/intrinsics/math.jl#L370)
    # This is adapted from Julia's mod implementation (also based on rem), but trying to
    # avoid branches (not sure it improves much).
    r = ifelse(iszero(r), copysign(r, L), r)  # replaces -0.0 -> +0.0
    # (2) Determine associated cell index.
    # Here unsafe_trunc(Int, ⋅) is used instead of floor(Int, ⋅) because it should be faster.
    # For non-negative values, both should give the same result.
    # The unsafe_trunc function generally uses a single intrinsic CPU instruction and never
    # throws errors. It can silently give a wrong result if the values are not representable
    # by an Int, but that will never be the case in practice here (since 0 ≤ x/rcut < L/rcut
    # and L/rcut is very small compared to typemax(Int) = 2^63 - 1).
    determine_cell_index_folded(x, rcut, N)
end

## ========================================================================================== ##

# Fallback implementation of atomicswap using CAS (compare-and-swap) loop, since Atomix.@atomicswap is not currently
# implemented (https://github.com/JuliaConcurrent/Atomix.jl/pull/62).
@inline function _atomicswap!(dst::AbstractArray{T}, val::T, idx, ::Val{order}) where {T <: Integer, order}
    # @inbounds old = Atomix.@atomicswap order dst[idx...] = val  # returns the old value // Not currently implemented in Atomix{CUDA,OpenCL,...}Ext
    success = false
    old = @inbounds dst[idx]
    while !success
        (; old, success) = @inbounds Atomix.@atomicreplace dst[idx] old => val
    end
    old
end

@kernel function set_elements_kernel!(
        get_coordinate::F, head_indices::PaddedArray, next_index, @Const(xp), rs_cell, Ls, folded::Val,
    ) where {F}
    n = @index(Global, Linear)
    x⃗ = @inline get_coordinate(@inbounds xp[n])  # usually get_coordinate === identity
    head_indices_data = parent(head_indices)   # full data associated to padded array
    nghosts = PaddedArrays.npad(head_indices)  # number of ghost cells per boundary (compile-time constant)
    if folded === Val(true)
        inds = map(determine_cell_index_folded, Tuple(x⃗), rs_cell, size(head_indices))
    else
        inds = map(determine_cell_index_gpu, Tuple(x⃗), rs_cell, Ls, size(head_indices))
    end
    I = CartesianIndex(inds .+ nghosts)  # shift by number of ghost cells, since we access raw data associated to padded array
    IndexType = eltype(head_indices_data)
    # @inbounds head_old = Atomix.@atomicswap :monotonic head_indices_data[I] = n  # returns the old value // Not currently implemented in AtomixCUDAExt
    head_old = _atomicswap!(head_indices_data, IndexType(n), I, Val(:monotonic))   # returns the old value
    @inbounds next_index[n] = head_old  # the old head now comes after the new element
    nothing
end

function _set_elements!(backend::GPU, get_coordinate::F, cl::PeriodicCellList, xp::AbstractVector, folded::Val) where {F}
    (; next_index, head_indices, Ls, rs_cell,) = cl
    head_indices_data = parent(head_indices)   # full data associated to padded array
    fill!(head_indices_data, EMPTY)
    Np = length(xp)
    resize!(next_index, Np)
    Base.require_one_based_indexing(xp)
    Base.require_one_based_indexing(next_index)
    groupsize = 256
    kernel! = set_elements_kernel!(backend, groupsize)
    kernel!(get_coordinate, head_indices, next_index, xp, rs_cell, Ls, folded; ndrange = size(xp))
    pad_periodic!(cl.head_indices)
    cl
end

## ========================================================================================== ##

@kernel function foreach_pair_kernel(
        f::F, @Const(head_indices::PaddedArray{M}), @Const(next_index), @Const(xp_dest),
        rs_cell, Ls,
        ::Val{M}, folded::Val,  # = number of subdivisions (equal to number of ghost cells)
    ) where {F <: Function, M}
    i = @index(Global, Linear)
    x⃗ = @inbounds xp_dest[i]
    if folded === Val(true)
        inds_central = map(determine_cell_index_folded, Tuple(x⃗), rs_cell, size(head_indices))
    else
        inds_central = map(determine_cell_index_gpu, Tuple(x⃗), rs_cell, Ls, size(head_indices))
    end
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    cell_indices = CartesianIndices(map(i -> (i - M):(i + M), Tuple(I₀)))
    @inbounds for I in cell_indices
        j = head_indices[I]
        while j != EMPTY
            @inline f(x⃗, i, j)
            j = next_index[j]
        end
    end
    nothing
end

@kernel function foreach_pair_batched_kernel(
        f::F, @Const(head_indices::PaddedArray{M}), @Const(next_index), @Const(xp_dest),
        rs_cell, Ls,
        ::Val{batch_size},
        ::Val{M}, folded::Val,  # = number of subdivisions (equal to number of ghost cells)
    ) where {F <: Function, batch_size, M}
    i = @index(Global, Linear)
    x⃗ = @inbounds xp_dest[i]
    IndexType = eltype(head_indices)
    inds = MVector{batch_size, IndexType}(undef)
    if folded === Val(true)
        inds_central = map(determine_cell_index_folded, Tuple(x⃗), rs_cell, size(head_indices))
    else
        inds_central = map(determine_cell_index_gpu, Tuple(x⃗), rs_cell, Ls, size(head_indices))
    end
    I₀ = CartesianIndex(inds_central)  # index of central cell (where x⃗ is located)
    cell_indices = CartesianIndices(map(i -> (i - M):(i + M), Tuple(I₀)))
    m = 0
    @inbounds for I in cell_indices
        j = head_indices[I]
        while j != EMPTY
            inds[m += 1] = j
            if m == batch_size
                @inline f(x⃗, i, Tuple(inds), batch_size)
                m = 0
            end
            j = next_index[j]
        end
    end
    if m > 0
        for l in (m + 1):batch_size
            # copy latest index, just to make sure that all returned indices are valid
            @inbounds inds[l] = inds[m]
        end
        @inline f(x⃗, i, Tuple(inds), m)
    end
    nothing
end

function _foreach_pair(backend::GPU, f::F, cl, xp_dest, batch_size, folded::Val, sort_points) where {F}
    (; head_indices, next_index, Ls, rs_cell,) = cl
    M = subdivisions(cl)
    groupsize = 256
    if batch_size === nothing
        kernel = foreach_pair_kernel(backend, groupsize)
        kernel(f, head_indices, next_index, xp_dest, rs_cell, Ls, Val(M), folded; ndrange = size(xp_dest))
    else
        kernel = foreach_pair_batched_kernel(backend, groupsize)
        kernel(f, head_indices, next_index, xp_dest, rs_cell, Ls, batch_size, Val(M), folded; ndrange = size(xp_dest))
    end
    nothing
end

