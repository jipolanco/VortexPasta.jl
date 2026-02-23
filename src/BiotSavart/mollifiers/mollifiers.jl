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

include("gaussian.jl")
