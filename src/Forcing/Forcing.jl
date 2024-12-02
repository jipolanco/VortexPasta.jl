"""
    Forcing

Defines methods for injecting energy onto a system of vortices.
"""
module Forcing

export AbstractForcing, NormalFluidForcing

using StaticArrays: SVector
using Random: Random, AbstractRNG
using LinearAlgebra: ⋅, ×

"""
    AbstractForcing

Abstract type representing a forcing method.
"""
abstract type AbstractForcing end

@doc raw"""
    NormalFluidForcing <: AbstractForcing

Abstract type representing a forcing method based on an imposed normal fluid flow which
modifies the vortex velocities via a mutual friction term.

This type of forcing defines an external velocity ``\bm{v}_{\text{f}}`` affecting vortex
motion, so that the total vortex velocity is

```math
\frac{\mathrm{d}\bm{s}}{\mathrm{d}t} = \bm{v}_{\text{s}} + \bm{v}_{\text{f}}.
```

Here ``\bm{v}_{\text{s}}`` is the self-induced velocity obtained by applying Biot–Savart's law.

The forcing velocity is of the form:

```math
\bm{v}_{\text{f}} = α \bm{s}' × (\bm{v}_{\text{n}} - \bm{v}_{\text{s}})
- α' \bm{s}' × \left[ \bm{s}' × (\bm{v}_{\text{n}} - \bm{v}_{\text{s}}) \right]
```

where ``\bm{s}'`` is the local unit tangent vector, and ``α`` and ``α'`` are non-dimensional
coefficients representing the intensity of Magnus and drag forces.

See also [`AbstractForcing`](@ref).
"""
abstract type NormalFluidForcing <: AbstractForcing end

include("fourier.jl")

end
