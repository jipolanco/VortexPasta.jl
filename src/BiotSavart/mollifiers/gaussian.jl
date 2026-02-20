export GaussianMollifier

"""
    GaussianMollifier{T <: AbstractFloat} <: AbstractMollifier
    GaussianMollifier(α::AbstractFloat)

Gaussian mollifier for Ewald summation.

This is the standard mollifier traditionally used in Ewald methods, associated to the
splitting `erf(αr) + erfc(αr) = 1`.

# Parameters

- `α::AbstractFloat`: splitting parameter (an inverse length scale).
"""
struct GaussianMollifier{T <: AbstractFloat} <: AbstractMollifier
    α::T
end

function Base.show(io::IO, g::GaussianMollifier)
    print(io, "GaussianMollifier(α = $(g.α))")
end

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

@inline exp_nosimd(::KA.Backend, x::Zero) = exp(x)    # = 1 (defined in Constants.jl)
@inline exp_simd(x::Zero) = exp(x)
@inline erf_nosimd(::KA.Backend, ::Zero) = Zero()
@inline erfc_nosimd(::KA.Backend, ::Zero) = One()
@inline erfc_simd(::Zero) = One()

@inline two_over_sqrt_pi(::SIMD.Vec{W, T}) where {W, T} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::T) where {T <: AbstractFloat} = 2 / sqrt(T(π))
@inline two_over_sqrt_pi(::Zero) = Zero()  # we don't really care about this value; it gets multiplied by Zero() anyway

@inline function weights_shortrange_simd(g::GaussianMollifier, r)
    (; α,) = g
    αr = α * r
    a = erfc_simd(αr)
    b = two_over_sqrt_pi(r) * αr * exp_simd(-(αr * αr))
    a, b
end

@inline function weights_shortrange_nosimd(backend::KA.Backend, g::GaussianMollifier, r)
    (; α,) = g
    αr = α * r
    a = erfc_nosimd(backend, αr)
    b = two_over_sqrt_pi(r) * αr * exp_nosimd(backend, -(αr * αr))
    a, b
end

# Same as above but replacing erfc -> erf.
@inline function weights_longrange_simd(g::GaussianMollifier, r)
    (; α,) = g
    αr = α * r
    a = erf_simd(αr)
    b = two_over_sqrt_pi(r) * αr * exp_simd(-(αr * αr))
    a, b
end

@inline function weights_longrange_nosimd(backend::KA.Backend, g::GaussianMollifier, r)
    (; α,) = g
    αr = α * r
    a = erf_nosimd(backend, αr)
    b = two_over_sqrt_pi(r) * αr * exp_nosimd(backend, -(αr * αr))
    a, b
end
