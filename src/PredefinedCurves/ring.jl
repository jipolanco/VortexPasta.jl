export Ring

@doc raw"""
    Ring(; r = Returns(0))

Describes a circular ring in 3D space.

By default, the ring has radius ``R = 1``, it is defined on the XY plane, and its centre is
located at `(0, 0, 0)`.

The circular ring can be optionally perturbed in the radial direction using a user-defined
function ``r(t)``, which must be periodic with period 1.
For example, one can pass `r(t) = 0.1 * cos(2πmt)` for a mode ``m`` sinusoidal
perturbation (associated to a wavelength ``λ = 2πR/m``) whose amplitude is 10% of the ring
radius.

Moreover, the geometry can be modified via the optional `transform` and `translate`
arguments of the [`define_curve`](@ref) function.

## Definition

```math
\begin{aligned}
    x(t) &= (1 + r(t)) cos(2πt) \\
    y(t) &= (1 + r(t)) sin(2πt) \\
    z(t) &= 0
\end{aligned}
```

for ``t ∈ [0, 1]``.

## Examples

Define a circular ring of radius 2:

```jldoctest PredefinedRing; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> p = Ring();

julia> S = define_curve(p; scale = 2);  # scale the curve by 2 (to get radius = 2)

julia> ts = range(0, 1; length = 16 + 1);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [2.0, 0.0, 0.0]
 [1.8477590650225735, 0.7653668647301796, 0.0]
 [1.4142135623730951, 1.414213562373095, 0.0]
 [0.7653668647301796, 1.8477590650225735, 0.0]
 [0.0, 2.0, 0.0]
 [-0.7653668647301796, 1.8477590650225735, 0.0]
 [-1.4142135623730951, 1.414213562373095, 0.0]
 [-1.8477590650225735, 0.7653668647301796, 0.0]
 [-2.0, 0.0, 0.0]
 [-1.8477590650225735, -0.7653668647301796, 0.0]
 [-1.4142135623730951, -1.414213562373095, 0.0]
 [-0.7653668647301796, -1.8477590650225735, 0.0]
 [0.0, -2.0, 0.0]
 [0.7653668647301796, -1.8477590650225735, 0.0]
 [1.4142135623730951, -1.414213562373095, 0.0]
 [1.8477590650225735, -0.7653668647301796, 0.0]
 [2.0, 0.0, 0.0]
```

Define a circular ring of radius 4 with a radial perturbation of 10%:

```jldoctest PredefinedRing; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> rfun(t) = 0.1 * cospi(4t);  # this is a mode m = 2 perturbation with 10% amplitude

julia> p = Ring(r = rfun);

julia> S = define_curve(p; scale = 4);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [4.4, 0.0, 0.0]
 [3.9568307230204227, 1.6389729494895986, 0.0]
 [2.8284271247461903, 2.82842712474619, 0.0]
 [1.4224945094311199, 3.4342055370698716, 0.0]
 [0.0, 3.6, 0.0]
 [-1.4224945094311199, 3.4342055370698716, 0.0]
 [-2.8284271247461903, 2.82842712474619, 0.0]
 [-3.9568307230204227, 1.6389729494895986, 0.0]
 [-4.4, 0.0, 0.0]
 [-3.9568307230204227, -1.6389729494895986, 0.0]
 [-2.8284271247461903, -2.82842712474619, 0.0]
 [-1.4224945094311199, -3.4342055370698716, 0.0]
 [0.0, -3.6, 0.0]
 [1.4224945094311199, -3.4342055370698716, 0.0]
 [2.8284271247461903, -2.82842712474619, 0.0]
 [3.9568307230204227, -1.6389729494895986, 0.0]
 [4.4, 0.0, 0.0]
```
"""
@kwdef struct Ring{FunR <: Function} <: ParametricCurve
    r :: FunR = Returns(0)
end

function _definition(p::Ring)
    (; r,) = p
    function S(t)
        u = 2 * t  # note: π * u ∈ [0, 2π]
        s, c = sincospi(u)
        R = 1 + r(t)
        SVector(R * c, R * s, 0)
    end
end
