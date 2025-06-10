"""
    helicity_spectrum(iter::VortexFilamentSolver; unfilter = true) -> (ks, Hk)
    helicity_spectrum(cache; unfilter = true) -> (ks, Hk)

Compute helicity spectrum associated to vortex filament state.

Returns a tuple of vectors `(ks, Hk)` where `ks` contains the probed wavenumbers and `Hk`
the helicity associated to each wavenumber.

See also [`helicity_spectrum!`](@ref) for a non-allocating variant and for more details.
"""
function helicity_spectrum end

"""
    helicity_spectrum!(Hk::AbstractVector, ks::AbstractVector, cache; unfilter = true)

Compute helicity spectrum associated to vortex filament state.

See [`energy_spectrum!`](@ref) for details on accepted `cache` arguments and on the optional
`unfilter` argument.
"""
function helicity_spectrum! end

function helicity_spectrum(cache; kws...)
    ks, Hk = init_spectrum(cache)
    helicity_spectrum!(Hk, ks, cache; kws...)
end

function helicity_spectrum!(
        Hk::AbstractVector, ks::AbstractVector, cache::LongRangeCache;
        unfilter = true,  # undo Ewald smoothing filter
    )
    (; state, ewald_op_d, ewald_prefactor,) = cache.common
    σ = BiotSavart.ewald_smoothing_scale(cache)
    from_smoothed_velocity = state.quantity == :velocity && state.smoothing_scale == σ  # smoothed velocity
    from_vorticity = state.quantity == :vorticity && state.smoothing_scale == 0  # unsmoothed vorticity
    γ² = ewald_prefactor^2  # = (Γ/V)^2
    if from_smoothed_velocity
        if unfilter
            _compute_spectrum!(Hk, ks, cache) do u⃗, k⃗, k², I
                # It's slightly faster to reuse values in ewald_op_d than to recompute exponentials...
                local w = @inbounds k² * ewald_op_d[I]
                local β = ifelse(
                    iszero(w),
                    one(γ²),   # set the factor to 1 if k² == 0
                    γ² / w^2,  # note: γ cancels out with prefactor already included in ewald_op_d
                )
                # @assert β ≈ exp(k² / (2 * cache.common.params.common.α^2))
                # @assert β ≈ exp(k² * σ^2)
                ω⃗ = im * (k⃗ × u⃗)
                H = u⃗ ⋅ ω⃗
                # @assert real(H) + imag(H) ≈ real(H)  # result is purely real
                # @assert H ≈ 2 * (k⃗ ⋅ (real(u⃗) × imag(u⃗)))  # alternative way of computing H(k⃗)
                2 * β * real(H)
            end
        else
            _compute_spectrum!(Hk, ks, cache) do u⃗, k⃗, k², I
                # Note: û ⋅ ω̂ is expected to be purely real.
                # Also, the factor 2 is because _compute_spectrum applies a factor 0.5 (for
                # energy), which is not present in the helicity.
                ω⃗ = im * (k⃗ × u⃗)
                2 * real(u⃗ ⋅ ω⃗)
            end
        end
    elseif from_vorticity
        _compute_spectrum!(Hk, ks, cache) do ω⃗, k⃗, k², I
            local u⃗ = im * (k⃗ × ω⃗) / k²
            2 * γ² * real(u⃗ × ω⃗)
        end
    else
        error(lazy"the current state of the long-range cache ($state) is currently not supported")
    end
    ks, Hk
end

helicity_spectrum!(Hk::AbstractVector, ks::AbstractVector, cache; kws...) =
    helicity_spectrum!(Hk, ks, get_long_range_cache(cache); kws...)
