"""
    PredefinedCurves

Module defining some commonly used parametric curves.
"""
module PredefinedCurves

export define_curve

abstract type ParametricCurve end

using LinearAlgebra: LinearAlgebra
using StaticArrays: SVector, SDiagonal

@doc raw"""
    define_curve(
        p::ParametricCurve;
        scale = 1, rotate = LinearAlgebra.I, translate = 0,
        orientation::Int = 1,
    ) -> Function

Return the definition of the parametric curve `p` as a function.

The returned function `S(t)` returns a coordinate ``\bm{x}`` for any given value of the scalar
parameter ``t ∈ [0, 1]``. In particular, closed curves satisfy `S(0) == S(1)`.

# Optional arguments

## Coordinate transformations

The original curve can be transformed by (1) scaling, (2) rotation and (3) translation
operations. Note that transformations are applied in that order.

- `scale`: scales the curve by a given factor. The argument can be a scalar value for
  isotropic scaling (same scaling in all directions) or a `Tuple` of values for anisotropic
  scaling (see below for some examples).

- `rotate`: in 3D, this should be a 3×3 orthogonal matrix describing pure rotation. For
  convenience, one can use the [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl)
  package for defining such rotations using different parametrisations (see below for some
  examples).

- `translate`: a scalar or a vector describing a translation operation. Note that
  translations are performed *after* scaling and rotation.

## Curve orientation

- `orientation`: allows to set the curve orientation. In particular, this determines the
  orientation of the tangent vector along the curve. Should be either `1` or `-1`.

# Extended help

## Examples

### Translated and scaled circular ring

Define a circular ring of radius ``R = 2`` centred at ``\bm{x}₀ = (0, 0, 1)`` and evaluate
its coordinates over equispaced points.

```jldoctest parametric_definition; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> S = define_curve(Ring(); translate = (0, 0, 1), scale = 2);

julia> ts = range(0, 1; length = 16 + 1)
0.0:0.0625:1.0

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [2.0, 0.0, 1.0]
 [1.8477590650225735, 0.7653668647301796, 1.0]
 [1.4142135623730951, 1.4142135623730951, 1.0]
 [0.7653668647301796, 1.8477590650225735, 1.0]
 [0.0, 2.0, 1.0]
 [-0.7653668647301796, 1.8477590650225735, 1.0]
 [-1.4142135623730951, 1.4142135623730951, 1.0]
 [-1.8477590650225735, 0.7653668647301796, 1.0]
 [-2.0, 0.0, 1.0]
 [-1.8477590650225735, -0.7653668647301796, 1.0]
 [-1.4142135623730951, -1.4142135623730951, 1.0]
 [-0.7653668647301796, -1.8477590650225735, 1.0]
 [0.0, -2.0, 1.0]
 [0.7653668647301796, -1.8477590650225735, 1.0]
 [1.4142135623730951, -1.4142135623730951, 1.0]
 [1.8477590650225735, -0.7653668647301796, 1.0]
 [2.0, 0.0, 1.0]
```

### Anisotropic scaling: ellipse from circular ring

If one wants an ellipse instead of a circle, one can simply apply an anisotropic scaling
transformation:

```jldoctest parametric_definition; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> using LinearAlgebra

julia> S = define_curve(Ring(); scale = (2, 1, 1));

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [2.0, 0.0, 0.0]
 [1.8477590650225735, 0.3826834323650898, 0.0]
 [1.4142135623730951, 0.7071067811865476, 0.0]
 [0.7653668647301796, 0.9238795325112867, 0.0]
 [0.0, 1.0, 0.0]
 [-0.7653668647301796, 0.9238795325112867, 0.0]
 [-1.4142135623730951, 0.7071067811865476, 0.0]
 [-1.8477590650225735, 0.3826834323650898, 0.0]
 [-2.0, 0.0, 0.0]
 [-1.8477590650225735, -0.3826834323650898, 0.0]
 [-1.4142135623730951, -0.7071067811865476, 0.0]
 [-0.7653668647301796, -0.9238795325112867, 0.0]
 [0.0, -1.0, 0.0]
 [0.7653668647301796, -0.9238795325112867, 0.0]
 [1.4142135623730951, -0.7071067811865476, 0.0]
 [1.8477590650225735, -0.3826834323650898, 0.0]
 [2.0, 0.0, 0.0]
```

### Rotated circular ring

If we wanted a circular ring defined on the YZ plane instead of the default XZ plane, we can
achieve this by rotating the original curve by 90° about the Y axis.

```jldoctest parametric_definition; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> using Rotations: RotY  # there's also RotX and RotZ

julia> rot = RotY(π / 2)  # rotation of 90° about the Y axis
3×3 RotY{Float64} with indices SOneTo(3)×SOneTo(3)(1.5708):
  6.12323e-17  0.0  1.0
  0.0          1.0  0.0
 -1.0          0.0  6.12323e-17

julia> rot = SMatrix(replace(x -> abs(x) < 1e-16 ? zero(x) : x, rot))  # remove spurious near-zero values
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  0.0  0.0  1.0
  0.0  1.0  0.0
 -1.0  0.0  0.0

julia> S = define_curve(Ring(); scale = 2, rotate = rot);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [0.0, 0.0, -2.0]
 [0.0, 0.7653668647301796, -1.8477590650225735]
 [0.0, 1.4142135623730951, -1.4142135623730951]
 [0.0, 1.8477590650225735, -0.7653668647301796]
 [0.0, 2.0, 0.0]
 [0.0, 1.8477590650225735, 0.7653668647301796]
 [0.0, 1.4142135623730951, 1.4142135623730951]
 [0.0, 0.7653668647301796, 1.8477590650225735]
 [0.0, 0.0, 2.0]
 [0.0, -0.7653668647301796, 1.8477590650225735]
 [0.0, -1.4142135623730951, 1.4142135623730951]
 [0.0, -1.8477590650225735, 0.7653668647301796]
 [0.0, -2.0, 0.0]
 [0.0, -1.8477590650225735, -0.7653668647301796]
 [0.0, -1.4142135623730951, -1.4142135623730951]
 [0.0, -0.7653668647301796, -1.8477590650225735]
 [0.0, 0.0, -2.0]
```

More generally, to rotate about an arbitrary axis `ê = [ex, ey, ez]` by an angle
`θ` (in radians), one can use `rot = AngleAxis(θ, ex, ey, ez)` from the Rotations.jl package.

### Random rotations

In addition, the Rotations.jl package allows to easily generate random and uniformly
distributed rotations:

```jldoctest parametric_definition; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> using Rotations: QuatRotation  # parametrise rotations using quaternions

julia> using Random: MersenneTwister

julia> rng = MersenneTwister(42);

julia> rot = rand(rng, QuatRotation)  # uniformly distributed random rotation
3×3 QuatRotation{Float64} with indices SOneTo(3)×SOneTo(3)(QuaternionF64(0.923391, -0.0602724, 0.307861, 0.22122)):
  0.712566  -0.445655   0.541886
  0.371433   0.894858   0.24752
 -0.59522    0.0249001  0.803177

julia> S = define_curve(Ring(); scale = 2, rotate = rot);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [1.4251329706813418, 0.7428666122301657, -1.190439083830269]
 [0.9755612628402096, 1.3712140745918031, -1.0807646293652284]
 [0.37746919621652514, 1.790806624183378, -0.806553557235094]
 [-0.27808913376434097, 1.9377650989455064, -0.4095520174421194]
 [-0.8913109140138611, 1.7897164032775457, 0.04980010440813476]
 [-1.3688386873583265, 1.3691996090301746, 0.501570611801321]
 [-1.637973179106087, 0.7402345861333226, 0.8769815402966743]
 [-1.6577411025987892, -0.00142444225909442, 1.1188799791393182]
 [-1.4251329706813418, -0.7428666122301657, 1.190439083830269]
 [-0.9755612628402096, -1.3712140745918031, 1.0807646293652284]
 [-0.37746919621652514, -1.790806624183378, 0.806553557235094]
 [0.27808913376434097, -1.9377650989455064, 0.4095520174421194]
 [0.8913109140138611, -1.7897164032775457, -0.04980010440813476]
 [1.3688386873583265, -1.3691996090301746, -0.501570611801321]
 [1.637973179106087, -0.7402345861333226, -0.8769815402966743]
 [1.6577411025987892, 0.00142444225909442, -1.1188799791393182]
 [1.4251329706813418, 0.7428666122301657, -1.190439083830269]
```
"""
function define_curve(
        p::ParametricCurve;
        translate = 0, scale = LinearAlgebra.I, rotate = LinearAlgebra.I,
        orientation::Int = 1,
    )
    S_base = _definition(p)
    scale_la = maybe_convert_scale(scale)
    A = rotate * scale_la  # linear transformation (we apply scaling first!)
    function S(t)
        τ = ifelse(orientation > 0, t, 1 - t)
        x⃗ = S_base(τ)
        translate .+ A * x⃗
    end
end

# Convert `scale` argument to something that can be used in linear algebra operations.
# This allows passing `scale = (sx, sy, sz)` to scale different directions differently.
# In that case the argument is converted to a diagonal matrix.
maybe_convert_scale(scale) = scale
maybe_convert_scale(scale::NTuple{3, Real}) = SDiagonal(scale)  # this is a static diagonal matrix

# Closed lines
include("ring.jl")
include("trefoil.jl")
include("lemniscate.jl")

# Periodic unclosed lines
include("periodic_line.jl")

end
