export PeriodicLine

@doc raw"""
    PeriodicLine(; x = Returns(0), y = Returns(0), r = Returns(0))

Describes an infinite (but unclosed) line which repeats itself periodically.

## Definition

By default, this defines a straight infinite line passing through the origin and oriented in
the Z direction.

The line can be perturbed along the X and/or Y directions using user-defined functions
``x(t)`` and ``y(t)`` which must be periodic with period ``T = 1``.
Alternative, one can pass a complex function ``r(t) = x(t) + im * y(t)``, which may be more
convenient in some cases.

The line is defined by:

```math
\begin{aligned}
    x(t) &= x(t) + ℜ[r(t)] \\
    y(t) &= y(t) + ℑ[r(t)] \\
    z(t) &= t - 1/2
\end{aligned}
```

for ``t ∈ [0, 1]``. Note that the line is "centred" at ``z = 0``, in the sense that
``z(1/2) = 0``.

Note that one can change the default period ``T = 1`` of the curve using the `scale`
argument of [`define_curve`](@ref). See below for some examples.

## Examples

Define a ``2π``-periodic curve with a sinusoidal perturbation along X:

```jldoctest PeriodicLine; filter = r"(\d*)\.(\d{14})\d+" => s"\1.\2***"
julia> xfun(t) = 0.1 * sinpi(2t)  # this function satisfies having period T = 1
xfun (generic function with 1 method)

julia> p = PeriodicLine(x = xfun);

julia> S = define_curve(p; scale = 2π);  # we scale the curve by 2π (this also scales the perturbation amplitude!)

julia> ts = range(0, 1; length = 16 + 1);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, -3.141592653589793]
 [0.24044709195373853, 0.0, -2.748893571891069]
 [0.4442882938158367, 0.0, -2.356194490192345]
 [0.5804906304278862, 0.0, -1.9634954084936207]
 [0.6283185307179586, 0.0, -1.5707963267948966]
 [0.5804906304278862, 0.0, -1.1780972450961724]
 [0.4442882938158367, 0.0, -0.7853981633974483]
 [0.24044709195373853, 0.0, -0.39269908169872414]
 [0.0, 0.0, 0.0]
 [-0.24044709195373853, 0.0, 0.39269908169872414]
 [-0.4442882938158367, 0.0, 0.7853981633974483]
 [-0.5804906304278862, 0.0, 1.1780972450961724]
 [-0.6283185307179586, 0.0, 1.5707963267948966]
 [-0.5804906304278862, 0.0, 1.9634954084936207]
 [-0.4442882938158367, 0.0, 2.356194490192345]
 [-0.24044709195373853, 0.0, 2.748893571891069]
 [0.0, 0.0, 3.141592653589793]
```

Note that the amplitude of the resulting curve is not the original ``0.1`` but instead
``0.1 × 2π = 0.6823…``. This is because passing `scale = 2π` scales the curve in all
directions, including the direction of the perturbation.

Instead, if one wanted to keep the original perturbation while extending the curve period
from ``1`` to ``2π``, one would need to apply scaling only along the Z direction:

```jldoctest PeriodicLine; filter = r"(\d*)\.(\d{14})\d+" => s"\1.\2***"
julia> using StaticArrays: SDiagonal

julia> S = define_curve(p; scale = SDiagonal(1, 1, 2π));

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, -3.141592653589793]
 [0.03826834323650898, 0.0, -2.748893571891069]
 [0.07071067811865477, 0.0, -2.356194490192345]
 [0.09238795325112868, 0.0, -1.9634954084936207]
 [0.1, 0.0, -1.5707963267948966]
 [0.09238795325112868, 0.0, -1.1780972450961724]
 [0.07071067811865477, 0.0, -0.7853981633974483]
 [0.03826834323650898, 0.0, -0.39269908169872414]
 [0.0, 0.0, 0.0]
 [-0.03826834323650898, 0.0, 0.39269908169872414]
 [-0.07071067811865477, 0.0, 0.7853981633974483]
 [-0.09238795325112868, 0.0, 1.1780972450961724]
 [-0.1, 0.0, 1.5707963267948966]
 [-0.09238795325112868, 0.0, 1.9634954084936207]
 [-0.07071067811865477, 0.0, 2.356194490192345]
 [-0.03826834323650898, 0.0, 2.748893571891069]
 [0.0, 0.0, 3.141592653589793]
```

Now the perturbation amplitude is ``0.1`` as we wanted.
"""
@kwdef struct PeriodicLine{FunX <: Function, FunY <: Function, FunR <: Function} <: ParametricCurve
    x :: FunX = Returns(0)
    y :: FunY = Returns(0)
    r :: FunR = Returns(0)
end

function _definition(p::PeriodicLine)
    (; x, y, r,) = p
    function (t)
        rt = r(t)
        SVector(
            x(t) + real(rt),
            y(t) + imag(rt),
            t - 1/2,
        )
    end
end
