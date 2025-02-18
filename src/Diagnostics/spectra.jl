export energy_spectrum, energy_spectrum!

"""
    energy_spectrum(iter::VortexFilamentSolver; unfilter = true) -> (ks, Ek)
    energy_spectrum(cache; unfilter = true) -> (ks, Ek)

Compute kinetic energy spectrum associated to vortex filament state.

Returns a tuple of vectors `(ks, Ek)` where `ks` contains the probed wavenumbers and `Ek`
the energy associated to each wavenumber.

See also [`energy_spectrum!`](@ref) for a non-allocating variant and for more details.
"""
function energy_spectrum end

"""
    energy_spectrum!(Ek::AbstractVector, ks::AbstractVector, cache; unfilter = true)

Compute kinetic energy spectrum associated to vortex filament state.

Here `cache` contains the results of long-range Biot–Savart computations. It can be either:

- a [`LongRangeCache`](@ref);
- a [`BiotSavartCache`](@ref) (which contains a `LongRangeCache`);
- a `VortexFilamentSolver` from the `Timestepping` module (which contains a `BiotSavartCache`).

The energy spectrum is computed from a recent Biot–Savart calculation using fast Ewald
summation. More precisely, it is computed from the long-range velocity field in Fourier
space. The [`LongRangeCache`](@ref) associated to the calculation is expected to currently
contain this field.

In its most usual state, a `LongRangeCache` contains the long-range velocity field in the
Ewald method, which is a Gaussian-filtered field (see e.g.
[`BiotSavart.to_smoothed_velocity!`](@ref)). By default this function undoes the Gaussian
filter, so that the returned kinetic energy spectrum is that of the unsmoothed velocity
(which is singular at vortex positions, so it presents a slow decay in wavenumber space).
One can pass `unfilter = false` to return the spectrum associated to the smoothed field.

The cache can also contain an unsmoothed vorticity field in Fourier space (the result
obtained right after performing a NUFFT from the filament locations, see
[`BiotSavart.compute_vorticity_fourier!`](@ref). In this case this
function does the right thing and also computes the spectrum of the associated (unsmoothed)
velocity field. Currently, the `unfilter` argument is ignored in this case.

The vectors `Ek` and `ks` are expected to have the same length. Moreover, the vector of
wavenumbers `ks` should satisfy `ks[begin] == 0` and have a constant step
`Δk = ks[i + 1] - ks[i]`. For convenience, the [`init_energy_spectrum`](@ref) function can
be used to create these vectors.

See also [`energy_spectrum`](@ref) for an allocating variant which doesn't need predefined
`Ek` and `ks` vectors.
"""
function energy_spectrum! end

get_long_range_cache(c::BiotSavartCache) = c.longrange

"""
    Diagnostics.init_energy_spectrum(cache) -> (ks, Ek)

Initialise fields for storing an energy spectrum.

Returns a wavenumber vector `ks` and an uninitialised energy spectrum `Ek` with the right
dimensions, which can be then passed to [`energy_spectrum!`](@ref).

The returned arrays are always on the CPU, even when the `cache` contains GPU data.

See [`energy_spectrum!`](@ref) for details on the `cache` argument.
"""
function init_energy_spectrum(cache::LongRangeCache)
    (; wavenumbers_d,) = cache.common
    kxs = adapt(Array, wavenumbers_d[1])  # make sure these are on the CPU
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)
    @assert with_hermitian_symmetry == (kxs[end] > 0)
    M = with_hermitian_symmetry ? length(kxs) : (length(kxs) + 1) ÷ 2
    @assert kxs[M] > 0
    Δk = kxs[begin + 1] - kxs[begin]
    ks = range(0, kxs[M]; step = Δk)
    Ek = similar(ks)  # this one is also on the CPU
    ks, Ek
end

init_energy_spectrum(c) = init_energy_spectrum(get_long_range_cache(c))

function energy_spectrum(cache; kws...)
    ks, Ek = init_energy_spectrum(cache)
    energy_spectrum!(Ek, ks, cache; kws...)
end

function energy_spectrum!(
        Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache;
        unfilter = true,  # undo Ewald smoothing filter
    )
    (; state, ewald_op_d, ewald_prefactor,) = cache.common
    from_smoothed_velocity = state.quantity == :velocity && state.smoothed
    from_vorticity = state.quantity == :vorticity && !state.smoothed
    γ² = ewald_prefactor^2  # = (Γ/V)^2
    if from_smoothed_velocity
        if unfilter
            energy_spectrum!(Ek, ks, cache) do u², k⃗, k², I
                # It's slightly faster to reuse values in ewald_op_d than to recompute exponentials...
                local w = @inbounds k² * ewald_op_d[I]
                β = ifelse(
                    iszero(w),
                    one(γ²),   # set the factor to 1 if k² == 0
                    γ² / w^2,  # note: γ cancels out with prefactor already included in ewald_op_d
                )
                # @assert β ≈ exp(k² / (2 * params.common.α^2))
                u² * β
            end
        else
            energy_spectrum!(Ek, ks, cache) do u², k⃗, k², I
                u²  # return unmodified coefficient
            end
        end
    elseif from_vorticity
        energy_spectrum!(Ek, ks, cache) do u², k⃗, k², I
            β = ifelse(iszero(k²), one(γ²), γ² / k²)
            β * u²
        end
    else
        error(lazy"the current state of the long-range cache ($state) is currently not supported")
    end
    ks, Ek
end

energy_spectrum!(Ek::AbstractVector, ks::AbstractVector, cache; kws...) =
    energy_spectrum!(Ek, ks, get_long_range_cache(cache); kws...)

energy_spectrum!(f::F, Ek::AbstractVector, ks::AbstractVector, cache; kws...) where {F <: Function} =
    energy_spectrum!(f, Ek, ks, get_long_range_cache(cache); kws...)

# This variant is for now internal and not documented.
# Note: `Ek` is generally an array on the CPU, while the spectrum is computed on the GPU if
# `cache` contains GPU data.
function energy_spectrum!(
        f::F, Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache,
    ) where {F <: Function}
    (; wavenumbers_d, uhat_d,) = cache.common

    backend = KA.get_backend(cache)  # CPU, GPU

    eachindex(ks) === eachindex(Ek) ||
        throw(DimensionMismatch("incompatible dimensions of vectors"))
    iszero(ks[begin]) || throw(ArgumentError("output wavenumbers should include k = 0"))
    Δk = ks[begin + 1] - ks[begin]  # we assume this is constant
    Δk_inv = 1 / Δk
    with_hermitian_symmetry = BiotSavart.has_real_to_complex(cache)

    if backend isa KA.CPU  # avoid assertion on the GPU, simply to avoid scalar indexing from the CPU (slow)
        kxs = wavenumbers_d[1]
        @assert with_hermitian_symmetry == (kxs[end] > 0)
    end

    # We pad the `ndrange` (the visited indices) so that all threads are active whitin each
    # workgroup. This makes it easier to implement the reduction operations needed to obtain
    # the energy spectra. However, it means that some indices `I` will be outside of the grid,
    # and this needs to be checked in the kernel.
    groupsize = ka_default_workgroupsize(backend, size(uhat_d)) :: NTuple
    ngroups = cld.(size(uhat_d), groupsize)  # number of workgroups in each direction
    ndrange = ngroups .* groupsize  # pad array dimensions
    Nk = length(Ek)
    uhat_comps = StructArrays.components(uhat_d)  # = (ux, uy, uz)
    kernel! = ka_generate_kernel(energy_spectrum_kernel!, backend, ndrange)

    # Use Bumper to allocate temporary CPU arrays (doesn't do anything for GPU arrays).
    buf = Bumper.default_buffer()
    @no_escape buf begin
        Ek_groups = if backend isa KA.CPU
            @alloc(eltype(Ek), Nk, ngroups...)
        else
            KA.allocate(backend, eltype(Ek), (Nk, ngroups...))
        end

        kernel!(f, Ek_groups, uhat_comps, wavenumbers_d, with_hermitian_symmetry, Δk_inv)  # execute kernel

        # Sum data from all workgroups, writing results onto Ek_d.
        Ek_d = if typeof(KA.get_backend(Ek)) === typeof(backend)
            fill!(Ek, 0)
            Ek  # output is already on the compute device
        else
            KA.zeros(backend, eltype(Ek), Nk)  # allocate array on the device
        end
        Base.mapreducedim!(identity, +, Ek_d, Ek_groups)  # sum values from all workgroups onto Ek_d

        if backend isa KA.GPU
            KA.unsafe_free!(Ek_groups)  # manually free GPU memory
        end
    end

    if Ek !== Ek_d
        copyto!(Ek, Ek_d)  # copy from device to host
        KA.unsafe_free!(Ek_d)
    end

    ks, Ek
end

# Computes one partial energy spectrum per workgroup.
# Here Ek_blocks is an array of dimensions (Nk, number_of_workgroups...).
# Note that `uhat` must be passed as a tuple of arrays (u, v, w) as using @Const with
# StructArrays currently fails for some reason...
# NOTE: this is probably not optimal for CPU threads.
@kernel function energy_spectrum_kernel!(
        f::F,
        Ek_blocks::AbstractArray, @Const(uhat::Tuple),
        @Const(wavenumbers::Tuple),
        @Const(with_hermitian_symmetry), @Const(Δk_inv),
    ) where {F}
    # Definitions which persist across @synchronize statements:
    @uniform begin
        T = eltype(Ek_blocks)
        Nk = size(Ek_blocks, 1)     # spectrum length
        Nt = prod(@groupsize())     # number of work items ("threads") per workgroup
        # Memory shared across work items (threads) in a single workgroup (block):
        E_sm = @localmem T Nt       # energy E(k) * dk for each visited k⃗ (one per work item)
        n_sm = @localmem UInt16 Nt  # bin index associated to each visited k⃗ (one per work item)
    end

    I = @index(Global, Cartesian)
    tid = @index(Local, Linear)  # thread id (in 1:Nt)

    # Step 1: compute pair (n, E) associated to this thread, where:
    #  - E = "energy" in the local k⃗
    #  - n = wavenumber bin corresponding to the local k⃗
    # Results are written to the shared memory arrays E_sm and n_sm (fast storage shared by the workgroup).
    if checkbounds(Bool, uhat[1], I)
        k⃗ = map((v, i) -> @inbounds(v[i]), wavenumbers, Tuple(I))
        kx = k⃗[1]

        # Find bin for current k⃗
        factor = ifelse(!with_hermitian_symmetry || iszero(kx), T(0.5), T(1.0))
        k² = sum(abs2, k⃗)
        knorm = sqrt(k²)
        n = 1 + unsafe_trunc(Int, knorm * Δk_inv + T(0.5))  # this implicitly assumes ks[begin] == 0

        # Compute energy at current k⃗ and fill local arrays
        if n ≤ Nk
            u² = zero(T)
            for d in eachindex(uhat)
                @inbounds u² += abs2(uhat[d][I])
            end
            v² = f(u², k⃗, k², I)  # possibly modifies the computed coefficient
            @inbounds E_sm[tid] = factor * v² * Δk_inv
            @inbounds n_sm[tid] = n
        else
            # If we're outside the wanted energy spectrum (k > kmax), we just add 0 energy to
            # some arbitrary bin.
            @inbounds E_sm[tid] = 0
            @inbounds n_sm[tid] = 1
        end
    else
        # This is if we're outside of the array indices.
        @inbounds E_sm[tid] = 0
        @inbounds n_sm[tid] = 1
    end

    @synchronize  # make sure E_sm and n_sm are synchronised across threads

    # Step 2: add results to global memory.
    group_idx = @index(Group, NTuple)
    @inbounds Ek_wg = view(Ek_blocks, 1:Nk, group_idx...)  # local part of the output array
    if tid == 1
        fill!(Ek_wg, 0)
        @inbounds for l ∈ 1:Nt
            E = E_sm[l]
            n = n_sm[l]
            Ek_wg[n] += E
        end
    end

    nothing
end
