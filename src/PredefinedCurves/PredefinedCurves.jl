"""
    PredefinedCurves

Module defining some commonly used parametric curves.
"""
module PredefinedCurves

export define_curve

abstract type ParametricCurve end

using LinearAlgebra: I  # default rotation
using StaticArrays: SVector

"""
    define_curve(
        p::ParametricCurve;
        orientation::Int = 1, transformation = 1, translation = 0,
    ) -> Function

Return the definition of the parametric curve `p` as a function.

The returned function `S(t)` returns a coordinate ``x⃗`` for any given value of the scalar
parameter ``t ∈ [0, 1]``. In particular, closed curves satisfy `S(0) == S(1)`.

## Optional arguments

- `orientation`: allows to set the curve orientation. In particular, this determines the
  orientation of the tangent vector along the curve. Should be either `1` or `-1`.

- `transformation`: *linear* coordinate transformation. In the most general case, this is a
  composition of a rotation and a scaling operation, described by a 3×3 transformation matrix.
  For pure isotropic scaling, this can be a scalar value. In the case of rotations, they can
  be conveniently defined using the
  [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl) package (see below for some
  examples).

- `translation`: a scalar or a vector describing a translation operation. Note that
  translations are performed *after* linear transformations such as rotations.

## Examples

### Translated and scaled circular ring

Define a circular ring of radius ``R = 2`` centred at ``x⃗₀ = (0, 0, 1)`` and evaluate
its coordinates over equispaced points.

```jldoctest parametric_definition
julia> S = define_curve(Ring(); translation = (0, 0, 1), transformation = 2);

julia> ts = range(0, 1; length = 16 + 1)
0.0:0.0625:1.0

julia> S.(ts)
17-element Vector{StaticArraysCore.SVector{3, Float64}}:
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

```jldoctest parametric_definition
julia> using LinearAlgebra, StaticArrays

julia> transformation = Diagonal(SVector(2.0, 1.0, 1.0))
3×3 Diagonal{Float64, SVector{3, Float64}} with indices SOneTo(3)×SOneTo(3):
 2.0   ⋅    ⋅
  ⋅   1.0   ⋅
  ⋅    ⋅   1.0

julia> S = define_curve(Ring(); transformation);

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

```jldoctest parametric_definition
julia> using Rotations: RotY  # there's also RotX and RotZ

julia> rot = RotY(π / 2)  # rotation of 90° about the Y axis
3×3 RotY{Float64} with indices SOneTo(3)×SOneTo(3)(1.5708):
  6.12323e-17  0.0  1.0
  0.0          1.0  0.0
 -1.0          0.0  6.12323e-17

julia> transformation = 2 * rot  # compose rotation with ×2 scaling (to get radius = 2)
3×3 SMatrix{3, 3, Float64, 9} with indices SOneTo(3)×SOneTo(3):
  1.22465e-16  0.0  2.0
  0.0          2.0  0.0
 -2.0          0.0  1.22465e-16

julia> S = define_curve(Ring(); transformation);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [1.2246467991473532e-16, 0.0, -2.0]
 [1.1314261122877003e-16, 0.7653668647301796, -1.8477590650225735]
 [8.659560562354934e-17, 1.4142135623730951, -1.4142135623730951]
 [4.6865204053262986e-17, 1.8477590650225735, -0.7653668647301796]
 [0.0, 2.0, 0.0]
 [-4.6865204053262986e-17, 1.8477590650225735, 0.7653668647301796]
 [-8.659560562354934e-17, 1.4142135623730951, 1.4142135623730951]
 [-1.1314261122877003e-16, 0.7653668647301796, 1.8477590650225735]
 [-1.2246467991473532e-16, 0.0, 2.0]
 [-1.1314261122877003e-16, -0.7653668647301796, 1.8477590650225735]
 [-8.659560562354934e-17, -1.4142135623730951, 1.4142135623730951]
 [-4.6865204053262986e-17, -1.8477590650225735, 0.7653668647301796]
 [0.0, -2.0, 0.0]
 [4.6865204053262986e-17, -1.8477590650225735, -0.7653668647301796]
 [8.659560562354934e-17, -1.4142135623730951, -1.4142135623730951]
 [1.1314261122877003e-16, -0.7653668647301796, -1.8477590650225735]
 [1.2246467991473532e-16, 0.0, -2.0]
```

More generally, to rotate about an arbitrary axis `ê = [ex, ey, ez]` by an angle
`θ` (in radians), one can use `rot = AngleAxis(θ, ex, ey, ez)` from the Rotations.jl package.

### Random rotations

In addition, the Rotations.jl package allows to easily generate random and uniformly
distributed rotations:

```jldoctest parametric_definition
julia> using Rotations: QuatRotation  # parametrise rotations using quaternions

julia> using Random: MersenneTwister

julia> rng = MersenneTwister(42);

julia> rot = rand(rng, QuatRotation)  # uniformly distributed random rotation
3×3 QuatRotation{Float64} with indices SOneTo(3)×SOneTo(3)(QuaternionF64(-0.719586, -0.575102, 0.0351433, -0.38758)):
 0.697094  -0.598216    0.395218
 0.517372   0.0380793  -0.854913
 0.496373   0.80043     0.336045

julia> S = define_curve(Ring(); transformation = 2 * rot);

julia> S.(ts)
17-element Vector{SVector{3, Float64}}:
 [1.394187993157782, 1.0347441397857304, 0.9927458914111738]
 [0.8302070621260415, 0.9851235420896431, 1.5297999348989406]
 [0.13983463173136992, 0.7855268152775545, 1.8339558059692869]
 [-0.5718263537403094, 0.4663407516577706, 1.8589085304315887]
 [-1.1964319604738316, 0.07615853598753197, 1.6008592821834717]
 [-1.6388916469079422, -0.32561812640796056, 1.0990937200484512]
 [-1.8318449366900935, -0.6778223807935073, 0.43000110234543676]
 [-1.7459164405768801, -0.9268343221784252, -0.3045552852199711]
 [-1.394187993157782, -1.0347441397857304, -0.9927458914111738]
 [-0.8302070621260415, -0.9851235420896431, -1.5297999348989406]
 [-0.13983463173136992, -0.7855268152775545, -1.8339558059692869]
 [0.5718263537403094, -0.4663407516577706, -1.8589085304315887]
 [1.1964319604738316, -0.07615853598753197, -1.6008592821834717]
 [1.6388916469079422, 0.32561812640796056, -1.0990937200484512]
 [1.8318449366900935, 0.6778223807935073, -0.43000110234543676]
 [1.7459164405768801, 0.9268343221784252, 0.3045552852199711]
 [1.394187993157782, 1.0347441397857304, 0.9927458914111738]
```
"""
function define_curve(
        p::ParametricCurve;
        orientation::Int = 1, translation = 0, transformation = 1,
    )
    S_base = _definition(p)
    function S(t)
        x⃗ = S_base(orientation * t)
        translation .+ transformation * x⃗
    end
end

include("ring.jl")

end
