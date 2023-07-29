export Lemniscate

@doc raw"""
    Lemniscate(; Az = 0,)

Describes a [lemniscate](https://en.wikipedia.org/wiki/Lemniscate) (or figure-eight curve)
in 3D space.

The specific definition used here corresponds to the [lemniscate of
Bernoulli](https://en.wikipedia.org/wiki/Lemniscate_of_Bernoulli) on the XY plane.
Under the normalised parametrisation ``x⃗(t)`` with ``t ∈ [0, 1]``, the curve crosses
the origin at ``t = 1/4`` and ``t = 3/4``, intersecting itself at that point.

One can perturb the curve in the third direction by passing a non-zero value of `Az` so that
the curve doesn't exactly intersect itself.

## Definition

```math
\begin{aligned}
    x(u) &= \frac{\cos u}{1 + \sin^2 u} \\
    y(u) &= \frac{\sin u \, \cos u}{1 + \sin^2 u} \\
    z(u) &= A_z \sin u
\end{aligned}
```

for ``u = 2πt ∈ [0, 2π]``.
"""
@kwdef struct Lemniscate{AmplitudeZ} <: ParametricCurve
    Az :: AmplitudeZ = 0
end

function lemniscate_bernoulli(t; Az)
    u = 2t  # in [0, 2]
    s, c = sincospi(u)
    denom = 1 + s * s
    SVector(
        c / denom,
        (s * c) / denom,
        Az * s,
    )
end

_definition(p::Lemniscate) = let Az = p.Az
    t -> lemniscate_bernoulli(t; Az)
end
