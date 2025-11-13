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
    clamp(1 + unsafe_trunc(Int, r / rcut), 1, N)  # make sure the index is in 1:N
end

@kernel function set_elements_kernel1!(
        get_coordinate::F, head_indices::PaddedArray, next_index, @Const(xp), rs_cell, Ls,
    ) where {F}
    n = @index(Global, Linear)
    x⃗ = @inline get_coordinate(@inbounds xp[n])  # usually get_coordinate === identity
    head_indices_data = parent(head_indices)   # full data associated to padded array
    nghosts = PaddedArrays.npad(head_indices)  # number of ghost cells per boundary (compile-time constant)
    inds = map(determine_cell_index_gpu, Tuple(x⃗), rs_cell, Ls, size(head_indices))
    I = CartesianIndex(inds .+ nghosts)  # shift by number of ghost cells, since we access raw data associated to padded array
    @inbounds head_old = Atomix.@atomicswap :monotonic head_indices_data[I] = n  # returns the old value
    @inbounds next_index[n] = head_old  # the old head now comes after the new element
    nothing
end

function _set_elements!(backend::GPU, get_coordinate::F, cl::PeriodicCellList, xp::AbstractVector) where {F}
    (; next_index, head_indices, Ls, rs_cell,) = cl
    head_indices_data = parent(head_indices)   # full data associated to padded array
    fill!(head_indices_data, EMPTY)
    Np = length(xp)
    resize!(next_index, Np)
    Base.require_one_based_indexing(xp)
    Base.require_one_based_indexing(next_index)
    groupsize = 256
    kernel! = set_elements_kernel1!(backend, groupsize)
    kernel!(get_coordinate, head_indices, next_index, xp, rs_cell, Ls; ndrange = size(xp))
    pad_periodic!(cl.head_indices)
    cl
end

