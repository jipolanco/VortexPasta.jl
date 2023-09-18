"""
    PredefinedCurves

Module defining some commonly used parametric curves.
"""
module PredefinedCurves

export define_curve

abstract type ParametricCurve end

using LinearAlgebra: LinearAlgebra
using StaticArrays: SVector

@doc raw"""
    define_curve(
        p::ParametricCurve;
        scale = 1, rotate = LinearAlgebra.I, translate = 0,
        orientation::Int = 1,
    ) -> Function

Return the definition of the parametric curve `p` as a function.

The returned function `S(t)` returns a coordinate ``x⃗`` for any given value of the scalar
parameter ``t ∈ [0, 1]``. In particular, closed curves satisfy `S(0) == S(1)`.

## Optional arguments

### Coordinate transformations

The original curve can be transformed by (1) scaling, (2) rotation and (3) translation
operations. Note that transformations are applied in that order.

- `scale`: scales the curve by a given factor. The argument can be a scalar value for
  isotropic scaling (same scaling in all directions), or a `Diagonal` matrix (3×3 in 3D) for
  anisotropic scaling. In the second case, one can use the `SDiagonal` type from the
  StaticArrays.jl package (see below for some examples).

- `rotate`: in 3D, this should be a 3×3 orthogonal matrix describing pure rotation. For
  convenience, one can use the [Rotations.jl](https://github.com/JuliaGeometry/Rotations.jl)
  package for defining such rotations using different parametrisations (see below for some
  examples).

- `translate`: a scalar or a vector describing a translation operation. Note that
  translations are performed *after* scaling and rotation.

### Curve orientation

- `orientation`: allows to set the curve orientation. In particular, this determines the
  orientation of the tangent vector along the curve. Should be either `1` or `-1`.

## Examples

### Translated and scaled circular ring

Define a circular ring of radius ``R = 2`` centred at ``x⃗₀ = (0, 0, 1)`` and evaluate
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
julia> using LinearAlgebra, StaticArrays

julia> scale = SDiagonal(2.0, 1.0, 1.0)
3×3 Diagonal{Float64, SVector{3, Float64}} with indices SOneTo(3)×SOneTo(3):
 2.0   ⋅    ⋅
  ⋅   1.0   ⋅
  ⋅    ⋅   1.0

julia> S = define_curve(Ring(); scale);

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
3×3 QuatRotation{Float64} with indices SOneTo(3)×SOneTo(3)(QuaternionF64(-0.719586, -0.575102, 0.0351433, -0.38758)):
 0.697094  -0.598216    0.395218
 0.517372   0.0380793  -0.854913
 0.496373   0.80043     0.336045

julia> S = define_curve(Ring(); scale = 2, rotate = rot);

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
        translate = 0, scale = LinearAlgebra.I, rotate = LinearAlgebra.I,
        orientation::Int = 1,
    )
    S_base = _definition(p)
    A = rotate * scale  # linear transformation (we apply scaling first!)
    function S(t)
        x⃗ = S_base(orientation * t)
        translate .+ A * x⃗
    end
end

# Closed lines
include("ring.jl")
include("trefoil.jl")
include("lemniscate.jl")

# Periodic unclosed lines
include("periodic_line.jl")

end
