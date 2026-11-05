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

# This inverts -∇²ψ = ω  <-->  k² ψ̂ = ω̂ (scalar version)
@kernel function fourier_inverse_laplacian_kernel_scalar!(us::AbstractArray, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = map((k, i) -> @inbounds(k[i]), ks, Tuple(I))
    k² = sum(abs2, k⃗)
    k²_inv = ifelse(iszero(k²), zero(k²), inv(k²))  # this implicitly sets the k⃗ = 0 mode to zero
    @inbounds us[I] *= k²_inv
    nothing
end

## ================================================================================ ##
## Adapt long-range backends to 2D point vortices

function _transform_to_fourier!(c::ExactSumCache{T, 2}, prefactor::Real) where {T}
    (; uhat, wavenumbers, pointdata,) = c.common
    (; charges,) = pointdata
    xp_src = get_source_points_longrange(pointdata)
    ω_hat = uhat.x  # vorticity is a scalar (uhat has 2 vector components)
    @assert size(ω_hat) == map(length, wavenumbers)
    fill!(ω_hat, zero(eltype(ω_hat)))
    inds = CartesianIndices(ω_hat)
    @inbounds for i in eachindex(xp_src, charges)
        X = xp_src[i]
        Q = prefactor * to_real(charges[i])
        @inbounds Threads.@threads for I in inds
            k⃗ = Vec2(map((ks, j) -> @inbounds(ks[j]), wavenumbers, Tuple(I)))
            ω_hat[I] += Q * cis(-k⃗ ⋅ X)
        end
    end
    # We zero out some "asymmetric" modes to ease the comparison with other implementations.
    _ensure_hermitian_symmetry!(c, ω_hat)
    c
end

# Case of a scalar field (usually the streamfunction in 2D point vortex)
function _interpolate_to_physical!(callback::F, output::AbstractVector{T}, c::ExactSumCache) where {F, T <: Real}
    (; wavenumbers, pointdata,) = c.common
    uhat = c.common.uhat.x  # scalar is contained in the first component of common.uhat
    xp_dest = get_destination_points_longrange(pointdata)
    @assert length(xp_dest) == length(output)
    kxs = first(wavenumbers)
    kx_lims = first(kxs), last(kxs)
    @assert kxs[2] > 0  # only positive half is included (Hermitian symmetry)
    @inbounds Threads.@threads for i in eachindex(xp_dest, output)
        X = xp_dest[i]
        q = zero(T)
        for I in CartesianIndices(uhat)
            k⃗ = SVector(map(getindex, wavenumbers, Tuple(I)))
            v_orig = uhat[I]::Complex{T}
            v = @inline callback(v_orig, Tuple(I))::Complex{T}
            z = cis(k⃗ ⋅ X)
            δq = if k⃗[1] in kx_lims
                # Note: the imaginary part will cancel out with -k⃗ (also computed)
                real(v) * real(z) - imag(v) * imag(z)  # = real(v * z)
            else
                # Apply Hermitian symmetry
                2 * (
                    real(v) * real(z) -
                    imag(v) * imag(z)
                ) # = real(v * z + conj(v) * conj(z))
            end
            q = q + δq
        end
        output[i] = q
    end
    nothing
end
