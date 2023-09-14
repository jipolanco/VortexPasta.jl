export energy_spectrum, energy_spectrum!

"""
    energy_spectrum(cache::LongRangeCache; unfilter = Val(true))      -> (ks, Ek)
    energy_spectrum(iter::VortexFilamentSolver; unfilter = Val(true)) -> (ks, Ek)

Compute kinetic energy spectrum associated to vortex filament state.

Returns a tuple of vectors `(ks, Ek)` where `ks` contains the probed wavenumbers and `Ek`
the energy associated to each wavenumber.

See also [`energy_spectrum!`](@ref) for a non-allocating variant and for more details.
"""
function energy_spectrum end

"""
    energy_spectrum!(Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache; unfilter = Val(true))
    energy_spectrum!(Ek::AbstractVector, ks::AbstractVector, iter::VortexFilamentSolver; unfilter = Val(true))

Compute kinetic energy spectrum associated to vortex filament state.

The energy spectrum is computed from a recent Biot–Savart calculation using fast Ewald
summation. More precisely, it is computed from the long-range velocity field in Fourier
space. The [`LongRangeCache`](@ref) associated to the calculation is expected to currently
contain this field.

Note that the long-range velocity field contained in a [`LongRangeCache`](@ref) is a
Gaussian-filtered field. By default this function undoes the Gaussian filter (unless one
passes `unfilter = Val(false)`), so that the spectrum displays a slow decay at large
wavenumbers (associated to the velocity being singular at vortex positions).

The vectors `Ek` and `ks` are expected to have the same length. Moreover, the vector of
wavenumbers `ks` should satisfy `ks[begin] == 0` and have a constant step
`Δk = ks[i + 1] - ks[i]`.

See also [`energy_spectrum`](@ref) for an allocating variant which doesn't need predefined
`Ek` and `ks` vectors.
"""
function energy_spectrum! end

function energy_spectrum(cache::LongRangeCache; kws...)
    (; wavenumbers,) = cache.common
    kxs = wavenumbers[1]
    with_hermitian_symmetry = kxs[end] > 0
    M = with_hermitian_symmetry ? length(kxs) : (length(kxs) + 1) ÷ 2
    @assert kxs[M] > 0
    Δk = kxs[begin + 1] - kxs[begin]
    ks_spectrum = range(0, kxs[M]; step = Δk)
    Ek = similar(ks_spectrum)
    energy_spectrum!(Ek, ks_spectrum, cache; kws...)
    ks_spectrum, Ek
end

function energy_spectrum!(
        Ek::AbstractVector, ks::AbstractVector, cache::LongRangeCache;
        unfilter::Val{UNFILTER} = Val(true),  # undo Ewald smoothing filter
    ) where {UNFILTER}
    (; state, wavenumbers, uhat, ewald_op, params,) = cache.common
    (; Γ, Ls,) = params.common
    eachindex(ks) === eachindex(Ek) ||
        throw(DimensionMismatch("incompatible dimensions of vectors"))
    state.quantity == :velocity && state.smoothed || error(
        lazy"the current state of the long-range cache ($state) is currently not supported"
    )
    iszero(ks[begin]) || throw(ArgumentError("output wavenumbers should include k = 0"))
    Δk = ks[begin + 1] - ks[begin]  # we assume this is constant
    kxs = wavenumbers[1]
    with_hermitian_symmetry = kxs[end] > 0
    fill!(Ek, 0)
    if UNFILTER
        ewald_prefactor = (Γ / prod(Ls))^2
    end
    @inbounds for I ∈ CartesianIndices(uhat)
        k⃗ = map(getindex, wavenumbers, Tuple(I))
        kx = k⃗[1]
        factor = (!with_hermitian_symmetry || kx == 0) ? 0.5 : 1.0
        k² = sum(abs2, k⃗)
        knorm = sqrt(k²)
        n = firstindex(Ek) + floor(Int, knorm / Δk)  # this implicitly assumes ks[begin] == 0
        n ≤ lastindex(Ek) || continue
        u² = sum(abs2, uhat[I])
        if UNFILTER
            # It's slightly faster to reuse values in ewald_op than to recompute
            # exponentials...
            β = ifelse(
                iszero(k²),
                one(ewald_prefactor),  # set the factor to 1 if k² == 0
                ewald_prefactor / (k² * ewald_op[I])^2,
            )
            # @assert β ≈ exp(k² / (2 * params.common.α^2))
            u² *= β
        end
        Ek[n] += factor * u²
    end
    Ek
end
