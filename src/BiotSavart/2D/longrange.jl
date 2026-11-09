function compute_vorticity_fourier!(cache::LongRangeCache{2})
    # Note: uhat has two components (it's a 2D vector field).
    # This function writes the Fourier-space vorticity to the first of these two components (uhat.x).
    (; uhat, state, vorticity_prefactor,) = cache.common
    transform_to_fourier!(cache, vorticity_prefactor)  # this writes the vorticity to uhat.x
    truncate_spherical!(cache)   # doesn't do anything if truncate_spherical = false (default)
    state.quantity = :vorticity
    state.smoothing_scale = 0
    uhat
end

function compute_streamfunction_fourier!(c::LongRangeCache{2})
    (; uhat, state,) = c.common
    wavenumbers = get_wavenumbers(c)
    from_vorticity = state.quantity === :vorticity && state.smoothing_scale == 0
    @assert from_vorticity
    ω_hat = uhat.x  # ω̂ is stored in the first component of uhat
    ka_backend = KA.get_backend(c)
    kernel = ka_generate_kernel(fourier_inverse_laplacian_kernel_scalar!, ka_backend, ω_hat)
    kernel(ω_hat, wavenumbers)
    state.quantity = :streamfunction
    uhat
end

function compute_velocity_fourier!(c::LongRangeCache{2})
    (; uhat, state,) = c.common
    wavenumbers = get_wavenumbers(c)
    from_vorticity = state.quantity === :vorticity && state.smoothing_scale == 0
    from_streamfunction = state.quantity === :streamfunction && state.smoothing_scale == 0
    @assert from_vorticity || from_streamfunction
    ka_backend = KA.get_backend(c)
    uhat_comps = StructArrays.components(uhat)
    if from_vorticity
        let kernel = ka_generate_kernel(fourier_biot_savart_2D_kernel!, ka_backend, uhat)
            kernel(uhat_comps, wavenumbers)
        end
    elseif from_streamfunction
        let kernel = ka_generate_kernel(fourier_curl_2D_kernel!, ka_backend, uhat)
            kernel(uhat_comps, wavenumbers)
        end
    end
    state.quantity = :velocity
    c
end

# This inverts -∇²ψ = ω  <-->  k² ψ̂ = ω̂ (scalar version)
@kernel function fourier_inverse_laplacian_kernel_scalar!(us::AbstractArray, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = map((k, i) -> @inbounds(k[i]), ks, Tuple(I))
    k² = sum(abs2, k⃗)
    k²_inv = ifelse(iszero(k²), zero(k²), inv(k²))  # this implicitly sets the k⃗ = 0 mode to zero
    @inbounds us[I] *= k²_inv
    nothing
end

# Applies 2D Biot-Savart operator in Fourier space:
# - at input, uhat[1] contains the vorticity field
# - at output, uhat[1:2] contains the velocity field
@kernel function fourier_biot_savart_2D_kernel!(uhat::NTuple{2}, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = map((k, i) -> @inbounds(k[i]), ks, Tuple(I))
    iω = @inbounds im * uhat[1][I]
    k² = sum(abs2, k⃗)
    k²_inv = ifelse(iszero(k²), zero(k²), inv(k²))  # this implicitly sets the k⃗ = 0 mode to zero
    @inbounds uhat[1][I] = +(k⃗[2] * k²_inv) * iω
    @inbounds uhat[2][I] = -(k⃗[1] * k²_inv) * iω
    nothing
end

# Applies 2D curl operator in Fourier space:
# - at input, uhat[1] contains the scalar field (typically streamfunction)
# - at output, uhat[1:2] contains the resulting vector field (typically velocity)
@kernel function fourier_curl_2D_kernel!(uhat::NTuple{2}, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = map((k, i) -> @inbounds(k[i]), ks, Tuple(I))
    iω = @inbounds im * uhat[1][I]
    @inbounds uhat[1][I] = +k⃗[2] * iω
    @inbounds uhat[2][I] = -k⃗[1] * iω
    nothing
end
