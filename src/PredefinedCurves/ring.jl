export Ring

"""
    Ring()

Describes a circular ring in 3D space.

By default, the ring has radius `1`, it is defined on the XY plane, and its centre is
located at `(0, 0, 0)`.

These default properties can be modified via the optional `transformation` and `translation`
arguments of the [`define_curve`](@ref) function.
"""
struct Ring <: ParametricCurve end

function _definition(::Ring)
    function S(t)
        u = 2 * t  # note: π * u ∈ [0, 2π]
        s, c = sincospi(u)
        SVector(c, s, 0)
    end
end
