export TrefoilKnot

@doc raw"""
    TrefoilKnot()

Describes a [trefoil knot](https://en.wikipedia.org/wiki/Trefoil_knot) in 3D space.

Before transformations, the trefoil knot is defined by:

```math
\begin{aligned}
    x &= \sin t + 2 \sin 2t \\
    y &= \cos t - 2 \cos 2t \\
    z &= -\sin 3t
\end{aligned}
```
"""
struct TrefoilKnot <: ParametricCurve end

function _definition(::TrefoilKnot)
    function S(t)
        u = 2t  # in [0, 2]
        s, c = sincospi(u)
        ss, cc = sincospi(2u)
        SVector(
            s + 2ss,
            c - 2cc,
            -3 * sinpi(3u),
        )
    end
end
