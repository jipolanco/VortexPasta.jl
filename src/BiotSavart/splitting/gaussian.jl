"""
    GaussianSplitting{T <: AbstractFloat} <: AbstractEwaldSplitting
    GaussianSplitting(; α::AbstractFloat)

Gaussian splitting kernel for Ewald summation.

This is the standard splitting kernel traditionally used in Ewald methods, associated to the
identity `erf(αr) + erfc(αr) = 1`.

# Parameters

- `α::AbstractFloat`: splitting parameter (an inverse length scale).
"""
@kwdef struct GaussianSplitting{T <: AbstractFloat} <: AbstractEwaldSplitting
    α::T
end

# This converts real values to the wanted precision.
convert_floats(::Type{T}, g::GaussianSplitting) where {T} = GaussianSplitting(convert(T, g.α))

function Base.show(io::IO, g::GaussianSplitting)
    print(io, "GaussianSplitting(α = $(g.α))")
end

accuracy_coefficient_shortrange(g::GaussianSplitting, rcut) = rcut / g.α
accuracy_coefficient_longrange(g::GaussianSplitting, kmax) = kmax / (2 * g.α)

# Evaluate splitting kernel in Fourier space.
# Note that this may be called from a GPU kernel.
# It doesn't need to be very performant since it's only done once when creating a BiotSavartCache.
function splitting_kernel_fourier(g::GaussianSplitting)
    (; α,) = g
    β = -1 / (4 * α * α)
    @inline(k² -> exp(β * k²))
end

# Factor in ⟨ψ⟩ = C * ⟨ω⟩ to be applied when the mean vorticity ⟨ω⟩ is nonzero.
# See background_vorticity_correction! for details.
background_vorticity_correction_factor(g::GaussianSplitting) = 1 / (4 * g.α^2)

include("SIMDFunctions/SIMDFunctions.jl")
using .SIMDFunctions: SIMDFunctions

@inline exp_simd(x::SIMD.Vec) = SIMDFunctions.exp(x)
@inline erf_simd(x::SIMD.Vec) = SIMDFunctions.erf(x)
@inline erfc_simd(x::SIMD.Vec) = SIMDFunctions.erfc(x)

# Note: even without explicit SIMD, calling SIMD-friendly implementations can enable
# automatic SIMD and thus noticeably improve performance.
@inline exp_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.exp(x)
@inline erf_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.erf(x)
@inline erfc_nosimd(::CPU, x::AbstractFloat) = SIMDFunctions.erfc(x)

# On GPU we call the functions from Base or SpecialFunctions, since these are usually
# overridden in each GPU implementation (CUDA, ...) and therefore should be fast.
@inline exp_nosimd(::GPU, x::AbstractFloat) = exp(x)
@inline erf_nosimd(::GPU, x::AbstractFloat) = SpecialFunctions.erf(x)
@inline erfc_nosimd(::GPU, x::AbstractFloat) = SpecialFunctions.erfc(x)

@inline two_over_sqrt_pi(::SIMD.Vec{W, T}) where {W, T} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))

@inline function weights_shortrange_simd(g::GaussianSplitting, r)
    (; α,) = g
    αr = α * r
    a = erfc_simd(αr)
    b = two_over_sqrt_pi(r) * αr * exp_simd(-(αr * αr))
    a, b
end

@inline function weights_shortrange_nosimd(backend::KA.Backend, g::GaussianSplitting, r)
    (; α,) = g
    αr = α * r
    a = erfc_nosimd(backend, αr)
    b = two_over_sqrt_pi(r) * αr * exp_nosimd(backend, -(αr * αr))
    a, b
end

# Same as above but replacing erfc -> erf.
# @inline function weights_longrange_simd(g::GaussianSplitting, r)
#     (; α,) = g
#     αr = α * r
#     a = erf_simd(αr)
#     b = two_over_sqrt_pi(r) * αr * exp_simd(-(αr * αr))
#     a, b
# end

@inline function weights_longrange_nosimd(backend::KA.Backend, g::GaussianSplitting, r)
    (; α,) = g
    αr = α * r
    a = erf_nosimd(backend, αr)
    b = two_over_sqrt_pi(r) * αr * exp_nosimd(backend, -(αr * αr))
    a, b
end
