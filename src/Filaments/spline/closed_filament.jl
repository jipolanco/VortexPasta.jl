export
    ClosedSplineFilament

@doc raw"""
    ClosedSplineFilament{T} <: ClosedFilament{T} <: AbstractFilament{T} <: AbstractVector{Vec3{T}}

Describes a closed curve (a loop) in 3D space using periodic cubic splines.

---

    ClosedSplineFilament(N::Integer, [T = Float64]; offset = zero(Vec3{T}))
    Filaments.init(ClosedFilament{T}, N::Integer, CubicSplineMethod(); offset = zero(Vec3{T}))

Allocate data for a closed spline filament with `N` discretisation points.

See also [`Filaments.init`](@ref).

# Examples

Initialise filament with set of discretisation points:

```jldoctest ClosedSplineFilament
julia> f = Filaments.init(ClosedFilament, 16, CubicSplineMethod());

julia> θs = range(-1, 1; length = 17)[1:16]
-1.0:0.125:0.875

julia> @. f = Vec3(cospi(θs), sinpi(θs), 0);

julia> f[4]
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.3826834323650898
 -0.9238795325112867
  0.0

julia> f[5] = (f[4] + 2 * f[6]) ./ 2
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  0.1913417161825449
 -1.38581929876693
  0.0

julia> update_coefficients!(f);
```

Note that [`update_coefficients!`](@ref) should be called whenever filament
coordinates are changed, before doing other operations such as estimating
derivatives.

Estimate derivatives at discretisation points:

```jldoctest ClosedSplineFilament
julia> f[4, Derivative(1)]
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  0.9090457394297016
 -0.727333461100651
  0.0

julia> f[4, Derivative(2)]
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  0.20911715113294102
 -2.09047051482799
  0.0
```

Estimate coordinates and derivatives in-between discretisation points:

```jldoctest ClosedSplineFilament
julia> f(4, 0.32)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.16753415613203387
 -1.1324592487590195
  0.0

julia> Ẋ, Ẍ = f(4, 0.32, Derivative(1)), f(4, 0.32, Derivative(2))
([0.8947546127964856, -0.9527970723463657, 0.0], [-0.3303413370703831, 0.17798009799460934, 0.0])

julia> X′, X″ = normalise_derivatives(Ẋ, Ẍ)
([0.6845546705034078, -0.7289615237390591, 0.0], [-0.05076295124082962, -0.047670575508846576, 0.0])
```

# Extended help

## Curve parametrisation

The parametrisation knots ``t_i`` are directly obtained from the interpolation point
positions.
A standard choice, which is used here, is for the knot increments to
approximate the arc length between two interpolation points:

```math
ℓ_{i} ≡ t_{i + 1} - t_{i} = |\bm{X}_{i + 1} - \bm{X}_i|,
```

which is a zero-th order approximation (and a lower bound) for the actual
arc length between points ``\bm{X}_i`` and ``\bm{X}_{i + 1}``.
"""
struct ClosedSplineFilament{
        T <: AbstractFloat,
        M,  # padding (fixed to 3?)
        Knots <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
    } <: ClosedFilament{T}

    # Spline knots: t_i = ∑_{j = 1}^{i - 1} ℓ_i
    ts :: Knots

    # Interpolation points X_i.
    Xs :: Points

    # Spline coefficients associated to the curve
    cs :: Points

    # Spline coefficients associated to first/second derivatives
    cderivs :: NTuple{2, Points}

    Xoffset :: Vec3{T}

    function ClosedSplineFilament(N::Integer, ::Type{T}; offset = zero(Vec3{T})) where {T}
        M = 3  # padding needed for cubic splines
        ts = PaddedVector{M}(Vector{T}(undef, N + 2M))
        Xs = similar(ts, Vec3{T})
        cs = similar(Xs)
        ċs = similar(Xs)
        c̈s = similar(Xs)
        cderivs = (ċs, c̈s)
        new{T, M, typeof(ts), typeof(Xs)}(ts, Xs, cs, cderivs, offset)
    end
end

init(::Type{ClosedFilament{T}}, N::Integer, ::CubicSplineMethod; kws...) where {T} =
    ClosedSplineFilament(N, T; kws...)

function Base.similar(f::ClosedSplineFilament, ::Type{T}, dims::Dims{1}) where {T <: Number}
    N, = dims
    ClosedSplineFilament(N, T; offset = f.Xoffset)
end

function Base.copyto!(v::ClosedSplineFilament, u::ClosedSplineFilament)
    copyto!(v.ts, u.ts)
    copyto!(v.Xs, u.Xs)
    copyto!(v.cs, u.cs)
    map(copyto!, v.cderivs, u.cderivs)
    v
end

discretisation_method(::ClosedSplineFilament) = CubicSplineMethod()
interpolation_method(::ClosedSplineFilament) = CubicSplineMethod()

function _update_coefficients_only!(f::ClosedSplineFilament)
    (; ts, Xs, cs, cderivs, Xoffset,) = f
    solve_cubic_spline_coefficients!(cs, ts, Xs; buf = cderivs[1], Xoffset,)
    spline_derivative!(cderivs[1], cs, ts, Val(4))
    spline_derivative!(cderivs[2], cderivs[1], ts, Val(3))
    f
end

(f::ClosedSplineFilament)(node::AtNode, ::Derivative{0} = Derivative(0)) = f[node.i]
(f::ClosedSplineFilament)(node::AtNode, der::Derivative) = f(node.i, 0.0, der)

# The second derivative at a node is simply equal to the corresponding spline coefficient.
(f::ClosedSplineFilament)(node::AtNode, ::Derivative{2}) = f.cderivs[2][node.i]

function (f::ClosedSplineFilament)(i::Int, ζ::Number, der::Derivative = Derivative(0))
    (; ts,) = f
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]  # convert ζ ∈ [0, 1] to spline parametrisation
    f(t, der; ileft = i)
end

# Here `t` is in the spline parametrisation.
function (f::ClosedSplineFilament)(
        t_in::Number, ::Derivative{n} = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    ) where {n}
    (; ts, cs, cderivs, Xoffset,) = f
    i, t = _find_knot_segment(ileft, knotlims(f), ts, t_in)
    coefs = (cs, cderivs...)[n + 1]
    ord = 4 - n
    y = evaluate_spline(coefs, ts, i, t, Val(ord))
    deperiodise_spline(y, Xoffset, ts, t_in, Val(n))  # only useful if Xoffset ≠ 0 ("infinite" / non-closed filaments)
end
