# Mollifiers (smoothing kernels) which can be used for performing Ewald summation.

"""
    AbstractMollifier

Abstract type representing a mollifier (smoothing kernel) which can be used for performing
Ewald summation.
"""
abstract type AbstractMollifier end

include("gaussian.jl")
