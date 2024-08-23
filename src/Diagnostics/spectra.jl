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

See [`energy_spectrum!`](@ref) for details on the `cache` argument.
"""
function init_energy_spectrum(cache::LongRangeCache)
    (; wavenumbers_d,) = cache.common
    kxs = wavenumbers_d[1]
    with_hermitian_symmetry = kxs[end] > 0
    M = with_hermitian_symmetry ? length(kxs) : (length(kxs) + 1) ÷ 2
    @assert kxs[M] > 0
    Δk = kxs[begin + 1] - kxs[begin]
    ks = range(0, kxs[M]; step = Δk)
    Ek = similar(ks)
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

# This variant is for now internal and not documented.
function energy_spectrum!(
        f::F, Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache,
    ) where {F <: Function}
    (; wavenumbers_d, uhat_d,) = cache.common
    eachindex(ks) === eachindex(Ek) ||
        throw(DimensionMismatch("incompatible dimensions of vectors"))
    iszero(ks[begin]) || throw(ArgumentError("output wavenumbers should include k = 0"))
    Δk = ks[begin + 1] - ks[begin]  # we assume this is constant
    Δk_inv = 1 / Δk
    kxs = wavenumbers_d[1]
    with_hermitian_symmetry = kxs[end] > 0
    fill!(Ek, 0)
    T = eltype(ks)
    @inbounds for I ∈ CartesianIndices(uhat_d)
        k⃗ = map((v, i) -> @inbounds(v[i]), wavenumbers_d, Tuple(I))
        kx = k⃗[1]
        factor = (!with_hermitian_symmetry || kx == 0) ? T(0.5) : T(1.0)
        k² = sum(abs2, k⃗)
        knorm = sqrt(k²)
        n = firstindex(Ek) + floor(Int, knorm * Δk_inv + T(0.5))  # this implicitly assumes ks[begin] == 0
        n ≤ lastindex(Ek) || continue
        u² = sum(abs2, uhat_d[I])
        v² = f(u², k⃗, k², I)  # possibly modifies the computed coefficient
        Ek[n] += factor * v² * Δk_inv
    end
    ks, Ek
end

energy_spectrum!(f::F, Ek::AbstractVector, ks::AbstractVector, cache; kws...) where {F <: Function} =
    energy_spectrum!(f, Ek, ks, get_long_range_cache(cache); kws...)
