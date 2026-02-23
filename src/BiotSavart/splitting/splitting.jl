# Smoothing kernels which can be used for performing Ewald summation.

export NoSplitting, GaussianSplitting

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

convert_floats(::Type{T}, g::NoSplitting) where {T} = g

Base.show(io::IO, ::NoSplitting) = print(io, "NoSplitting()")

accuracy_coefficient_shortrange(::NoSplitting, rcut) = rcut === Infinity() ? Infinity() : Zero()
accuracy_coefficient_longrange(::NoSplitting, kmax) = Infinity()  # there's no long-range

@inline weights_shortrange_simd(::NoSplitting, r) = one(r), zero(r)
@inline weights_shortrange_nosimd(::KA.Backend, ::NoSplitting, r) = one(r), zero(r)

# There's no long-range when splitting is disabled.
# @inline weights_longrange_simd(g::NoSplitting, r) = zero(r), zero(r)
# @inline weights_longrange_nosimd(::KA.Backend, g::NoSplitting, r) = zero(r), zero(r)

include("gaussian.jl")
