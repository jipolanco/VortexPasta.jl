# Smoothing kernels which can be used for performing Ewald summation.

export NoSplitting, GaussianSplitting, KaiserBesselSplitting

"""
    AbstractEwaldSplitting

Abstract type representing a splitting kernel used to perform Ewald summation.
"""
abstract type AbstractEwaldSplitting end

"""
    NoSplitting <: AbstractEwaldSplitting
    NoSplitting()

Represents the absence of an Ewald splitting kernel.

This may be used to completely disable Ewald summation and to impose open (non-periodic)
boundary conditions.

When choosing `NoSplitting()`, a naive (slow) method is used to compute Biot--Savart
interactions between filaments. Therefore, it is not recommended to use this for simulation
of dense vortex systems.
"""
struct NoSplitting <: AbstractEwaldSplitting end

periods(g::NoSplitting) = (Infinity(), Infinity(), Infinity())
cutoff_distance(g::NoSplitting) = Infinity()
fourier_grid_size(g::NoSplitting) = (0, 0, 0)

convert_floats(::Type{T}, g::NoSplitting) where {T} = g

Base.show(io::IO, ::NoSplitting) = print(io, "NoSplitting()")
Base.summary(io::IO, g::NoSplitting) = show(io, g)

@inline weights_shortrange_simd(::NoSplitting, r) = one(r), zero(r)
@inline weights_shortrange_nosimd(::KA.Backend, ::NoSplitting, r) = one(r), zero(r)

# There's no long-range when splitting is disabled.
# @inline weights_longrange_simd(g::NoSplitting, r) = zero(r), zero(r)
# @inline weights_longrange_nosimd(::KA.Backend, g::NoSplitting, r) = zero(r), zero(r)

include("chebyshev.jl")  # for function approximation (integral of KB kernel, ...)
using .ChebyshevApproximations: ChebyshevApproximations, ChebyshevSeries

include("gaussian.jl")
include("kaiser_bessel.jl")
