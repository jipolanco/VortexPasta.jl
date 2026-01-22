export energy_spectrum, energy_spectrum!

"""
    energy_spectrum(iter::VortexFilamentSolver) -> (ks, Ek)
    energy_spectrum(cache) -> (ks, Ek)

Compute kinetic energy spectrum associated to vortex filament state.

Returns a tuple of vectors `(ks, Ek)` where `ks` contains the probed wavenumbers and `Ek`
the energy associated to each wavenumber.

See also [`energy_spectrum!`](@ref) for a non-allocating variant and for more details.
"""
function energy_spectrum end

"""
    energy_spectrum!(Ek::AbstractVector, ks::AbstractVector, cache)

Compute kinetic energy spectrum associated to vortex filament state.

Here `cache` contains the results of long-range Biot–Savart computations. It can be either:

- a [`LongRangeCache`](@ref);
- a [`BiotSavartCache`](@ref) (which contains a `LongRangeCache`);
- a `VortexFilamentSolver` from the `Timestepping` module (which contains a `BiotSavartCache`).

The energy spectrum is computed from a recent Biot–Savart calculation using fast Ewald
summation. More precisely, it is computed from the (Fourier-truncated) velocity field in
Fourier space. The [`LongRangeCache`](@ref) associated to the calculation is expected to
currently contain this field.

The cache can also contain the vorticity field in Fourier space (the result
obtained right after performing a NUFFT from the filament locations, see
[`BiotSavart.compute_vorticity_fourier!`](@ref). In this case this function does the right
thing and also computes the spectrum of the associated velocity field.

The vectors `Ek` and `ks` are expected to have the same length. Moreover, the vector of
wavenumbers `ks` should satisfy `ks[begin] == 0` and have a constant step
`Δk = ks[i + 1] - ks[i]`. For convenience, the [`init_spectrum`](@ref) function can
be used to create these vectors.

See also [`energy_spectrum`](@ref) for an allocating variant which doesn't need predefined
`Ek` and `ks` vectors.
"""
function energy_spectrum! end

get_long_range_cache(c::BiotSavartCache) = c.longrange

function init_spectrum_wavenumbers(cache::LongRangeCache)
    (; wavenumbers,) = BiotSavart.get_longrange_field_fourier(cache)
    kxs = adapt(Array, wavenumbers[1])  # make sure these are on the CPU
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)
    @assert with_hermitian_symmetry == (kxs[end] > 0)
    M = with_hermitian_symmetry ? length(kxs) : (length(kxs) + 1) ÷ 2
    @assert kxs[M] > 0
    @assert kxs[begin] == 0
    Δk = kxs[begin + 1] - kxs[begin]
    range(0, kxs[M]; step = Δk)
end

"""
    Diagnostics.init_spectrum(iter::VortexFilamentSolver) -> (ks, Ek)
    Diagnostics.init_spectrum(cache) -> (ks, Ek)

Initialise fields for storing an energy or helicity spectrum.

Returns a wavenumber vector `ks` and an uninitialised spectrum `Ek` with the right
dimensions, which can be then passed to [`energy_spectrum!`](@ref) or [`helicity_spectrum!](@ref).

The returned arrays are always on the CPU, even when the `cache` contains GPU data.

See [`energy_spectrum!`](@ref) for details on the `cache` argument.
"""
function init_spectrum(cache::LongRangeCache)
    ks = init_spectrum_wavenumbers(cache)
    Ek = similar(ks)  # this one is also on the CPU
    ks, Ek
end

init_spectrum(c) = init_spectrum(get_long_range_cache(c))

init_energy_spectrum(c) = init_spectrum(c)  # for backwards compatibility

function energy_spectrum(cache; kws...)
    ks, Ek = init_spectrum(cache)
    energy_spectrum!(Ek, ks, cache; kws...)
end

function energy_spectrum!(
        Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache,
    )
    (; state,) = cache.common
    from_velocity = state.quantity == :velocity && state.smoothing_scale == 0
    from_vorticity = state.quantity == :vorticity && state.smoothing_scale == 0
    if from_velocity
        _compute_spectrum!(Ek, ks, cache) do u⃗, k⃗, k², I
            sum(abs2, u⃗)  # return unmodified coefficient
        end
    elseif from_vorticity
        _compute_spectrum!(Ek, ks, cache) do u⃗, k⃗, k², I
            local u² = sum(abs2, u⃗)
            local β = ifelse(iszero(k²), one(k²), 1 / k²)
            β * u²
        end
    else
        error(lazy"the current state of the long-range cache ($state) is currently not supported")
    end
    ks, Ek
end

energy_spectrum!(Ek::AbstractVector, ks::AbstractVector, cache; kws...) =
    energy_spectrum!(Ek, ks, get_long_range_cache(cache); kws...)

# This variant is for now internal and not documented.
# Note: `Ek` is generally an array on the CPU, while the spectrum is computed on the GPU if
# `cache` contains GPU data.
function _compute_spectrum!(
        f::F, Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache,
    ) where {F <: Function}
    (; field, wavenumbers,) = BiotSavart.get_longrange_field_fourier(cache)
    uhat_comps = field::NTuple  # = (ux, uy, uz)

    backend = KA.get_backend(cache)  # CPU, GPU

    eachindex(ks) === eachindex(Ek) ||
        throw(DimensionMismatch("incompatible dimensions of vectors"))
    iszero(ks[begin]) || throw(ArgumentError("output wavenumbers should include k = 0"))
    backend = KA.get_backend(cache)  # CPU, GPU
    _compute_spectrum_impl!(f, backend, Ek, ks, cache)
end

function _compute_spectrum_impl!(f::F, backend::KA.CPU, Ek, ks, cache) where {F <: Function}
    (; field, wavenumbers,) = BiotSavart.get_longrange_field_fourier(cache)
    uhat = field::NTuple  # = (ux, uy, uz)
    Δk = ks[begin + 1] - ks[begin]  # we assume this is constant
    Δk_inv = 1 / Δk
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)
    kxs = wavenumbers[1]
    @assert with_hermitian_symmetry == (kxs[end] > 0)

    Nk = length(Ek)

    # Use Bumper to allocate temporary CPU arrays.
    buf = Bumper.default_buffer()
    @no_escape buf begin
        nchunks = Threads.nthreads()
        T = eltype(Ek)
        Ek_threads = @alloc(T, Nk, nchunks)
        Base.require_one_based_indexing(Ek_threads)
        inds_all = CartesianIndices(uhat[1])
        @sync for j in 1:nchunks
            Threads.@spawn begin
                Ek_local = @view Ek_threads[:, j]
                fill!(Ek_local, zero(T))
                Nall = length(inds_all)
                istart = ((j - 1) * Nall) ÷ nchunks + 1
                iend = (j * Nall) ÷ nchunks
                inds = inds_all[istart:iend]
                @inbounds for I in inds
                    k⃗ = Vec3(map((v, i) -> @inbounds(v[i]), wavenumbers, Tuple(I)))
                    kx = k⃗[1]::T

                    # Find bin for current k⃗
                    factor = ifelse(!with_hermitian_symmetry || iszero(kx), T(0.5), T(1.0))
                    k² = sum(abs2, k⃗)::T
                    knorm = sqrt(k²)
                    n = 1 + unsafe_trunc(Int, knorm * Δk_inv + T(0.5))  # this implicitly assumes ks[begin] == 0

                    if n ≤ Nk
                        u⃗ = Vec3(
                            ntuple(Val(3)) do d
                                @inline
                                @inbounds(uhat[d][I])::Complex{T}
                            end
                        )
                        v² = f(u⃗, k⃗, k², I)  # possibly modifies the input coefficients
                        Ek_local[n] += factor * v² * Δk_inv
                    end
                end
            end
        end
        fill!(Ek, zero(eltype(Ek)))
        Base.mapreducedim!(identity, +, Ek, Ek_threads)  # sum values from all workgroups onto Ek
    end

    ks, Ek
end

function _compute_spectrum_impl!(f::F, backend::KA.GPU, Ek, ks, cache) where {F <: Function}
    (; field, wavenumbers,) = BiotSavart.get_longrange_field_fourier(cache)
    uhat_comps = field::NTuple  # = (ux, uy, uz)
    Δk = ks[begin + 1] - ks[begin]  # we assume this is constant
    Δk_inv = 1 / Δk
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)
    kxs = wavenumbers[1]
    @assert with_hermitian_symmetry == (kxs[end] > 0)

    # We pad the `ndrange` (the visited indices) so that all threads are active whitin each
    # workgroup. This makes it easier to implement the reduction operations needed to obtain
    # the energy spectra. However, it means that some indices `I` will be outside of the grid,
    # and this needs to be checked in the kernel.
    dims = size(uhat_comps[1])
    groupsize = (8, 8, 4)
    ngroups = cld.(dims, groupsize)  # number of workgroups in each direction
    ndrange = ngroups .* groupsize   # pad array dimensions
    Nk = length(Ek)
    kernel! = compute_spectrum_kernel_gpu!(backend, groupsize, ndrange)

    Ek_groups = KA.allocate(backend, eltype(Ek), (Nk, ngroups...))
    fill!(Ek_groups, zero(eltype(Ek_groups)))
    kernel!(f, Ek_groups, uhat_comps, wavenumbers, with_hermitian_symmetry, Δk_inv)  # execute kernel
    Ek_to_reduce = reshape(Ek_groups, (Nk, :))
    Ek_d = AK.reduce(+, Ek_to_reduce; init = zero(eltype(Ek_to_reduce)), dims = 2)  # sum values from all workgroups onto Ek_d

    copyto!(Ek, Ek_d)  # copy from device to host (or device to device, if Ek is on the same device as Ek_d)
    KA.synchronize(backend)  # just in case, to avoid freeing memory still being used
    KA.unsafe_free!(Ek_d)  # manually free GPU memory
    KA.unsafe_free!(Ek_groups)

    ks, Ek
end

# Computes one partial energy spectrum per workgroup.
# Here Ek_blocks is an array of dimensions (Nk, number_of_workgroups...).
# Note that `uhat` must be passed as a tuple of arrays (u, v, w).
@kernel unsafe_indices=true cpu=false function compute_spectrum_kernel_gpu!(
        f::F,
        Ek_blocks::AbstractArray, @Const(uhat::Tuple),
        @Const(wavenumbers::Tuple),
        @Const(with_hermitian_symmetry), @Const(Δk_inv),
    ) where {F}
    T = eltype(Ek_blocks)
    Nk = size(Ek_blocks, 1)     # spectrum length

    # We use unsafe_indices=true since we explicitly check that we're inbounds below.
    # This means that we should avoid explicitly asking for global indices.
    idx_local = @index(Local, NTuple)
    idx_group = @index(Group, NTuple)
    group_size = @groupsize()
    I = CartesianIndex(
        (@. idx_local + (idx_group - 1) * group_size)
    )

    tid = @index(Local, Linear)  # thread id (in 1:Nt)

    # Step 1: compute pair (n, E) associated to this thread, where:
    #  - E = "energy" in the local k⃗
    #  - n = wavenumber bin corresponding to the local k⃗
    # Results are written to the shared memory arrays E_sm and n_sm (fast storage shared by the workgroup).
    if checkbounds(Bool, uhat[1], I)
        k⃗ = Vec3(map((v, i) -> @inbounds(v[i]), wavenumbers, Tuple(I)))
        kx = k⃗[1]::T

        # Find bin for current k⃗
        factor = ifelse(!with_hermitian_symmetry || iszero(kx), T(0.5), T(1.0))
        k² = sum(abs2, k⃗)::T
        knorm = sqrt(k²)
        n = 1 + unsafe_trunc(Int, knorm * Δk_inv + T(0.5))  # this implicitly assumes ks[begin] == 0

        # Compute energy at current k⃗ and fill local arrays
        if n ≤ Nk
            u⃗ = Vec3(
                ntuple(Val(3)) do d
                    @inline
                    @inbounds(uhat[d][I])::Complex{T}
                end
            )
            v² = @inline f(u⃗, k⃗, k², I)  # possibly modifies the input coefficients
            E = factor * v² * Δk_inv
            @inbounds Atomix.@atomic Ek_blocks[n, idx_group...] += E
        end
    end

    nothing
end
