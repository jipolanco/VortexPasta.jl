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

```jldoctest ClosedSplineFilament; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> f = Filaments.init(ClosedFilament, 16, CubicSplineMethod());

julia> θs = range(-1, 1; length = 17)[1:16]
-1.0:0.125:0.875

julia> @. f = Vec3(cospi(θs), sinpi(θs), 0);

julia> f[4]
3-element SVector{3, Float64} with indices SOneTo(3):
 -0.3826834323650898
 -0.9238795325112867
  0.0

julia> f[5] = (f[4] + 2 * f[6]) ./ 2
3-element SVector{3, Float64} with indices SOneTo(3):
  0.1913417161825449
 -1.38581929876693
  0.0

julia> update_coefficients!(f);
```

Note that [`update_coefficients!`](@ref) should be called whenever filament
coordinates are changed, before doing other operations such as estimating
derivatives.

Estimate derivatives at discretisation points:

```jldoctest ClosedSplineFilament; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> f[4, Derivative(1)]
3-element SVector{3, Float64} with indices SOneTo(3):
  0.9090457394297018
 -0.7273334611006509
  0.0

julia> f[4, Derivative(2)]
3-element SVector{3, Float64} with indices SOneTo(3):
  0.20911715113294102
 -2.09047051482799
  0.0
```

Estimate coordinates and derivatives in-between discretisation points:

```jldoctest ClosedSplineFilament; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
julia> f(4, 0.32)
3-element SVector{3, Float64} with indices SOneTo(3):
 -0.16753415613203387
 -1.1324592487590195
  0.0

julia> Ẋ, Ẍ = f(4, 0.32, Derivative(1)), f(4, 0.32, Derivative(2))
([0.8947546127964856, -0.9527970723463657, 0.0], [-0.3303413370703831, 0.17798009799460934, 0.0])

julia> X′, X″ = f(4, 0.32, UnitTangent()), f(4, 0.32, CurvatureVector())
([0.6845546705034081, -0.7289615237390588, 0.0], [-0.050762951240829336, -0.047670575508846375, 0.0])
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
        M,  # padding (fixed to `npad(CubicSplineMethod)` == 3)
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
end

function ClosedSplineFilament(
        Xs::PaddedVector{M, Vec3{T}}; offset = zero(Vec3{T})
    ) where {M, T}
    @assert M == npad(discretisation_method(ClosedSplineFilament))
    ts = similar(Xs, T)
    cs = similar(Xs)
    cderivs = (similar(Xs), similar(Xs))
    Xoffset = convert(Vec3{T}, offset)
    ClosedSplineFilament(ts, Xs, cs, cderivs, Xoffset)
end

_init_closed_filament(Xs, ::CubicSplineMethod; kws...) = ClosedSplineFilament(Xs; kws...)

function change_offset(f::ClosedSplineFilament{T}, offset::Vec3) where {T}
    Xoffset = convert(Vec3{T}, offset)
    ClosedSplineFilament(f.ts, f.Xs, f.cs, f.cderivs, Xoffset)
end

allvectors(f::ClosedSplineFilament) = (f.ts, f.Xs, f.cs, f.cderivs...)

function Base.similar(f::ClosedSplineFilament, ::Type{T}, dims::Dims{1}) where {T <: Number}
    Xs = similar(nodes(f), Vec3{T}, dims)
    ClosedSplineFilament(Xs; offset = f.Xoffset)
end

discretisation_method(::Type{<:ClosedSplineFilament}) = CubicSplineMethod()
discretisation_method(::ClosedSplineFilament) = CubicSplineMethod()
interpolation_method(::ClosedSplineFilament) = CubicSplineMethod()

function _update_coefficients_only!(f::ClosedSplineFilament; only_derivatives = false)
    (; ts, Xs, cs, cderivs, Xoffset,) = f
    if !only_derivatives
        solve_cubic_spline_coefficients!(cs, ts, Xs; Xoffset,)
    end
    spline_derivative!(cderivs[1], cs, ts, Val(4))
    spline_derivative!(cderivs[2], cderivs[1], ts, Val(3))
    f
end

(f::ClosedSplineFilament)(node::AtNode, ::Derivative{0} = Derivative(0)) = f[node.i]

# The first derivative is a quadratic spline (order k = 3), and the formula to evaluate it
# on a knot is quite simple.
# This should give the same result as f(node.i, 0.0, Derivative(1)), but slightly faster.
function (f::ClosedSplineFilament)(node::AtNode, ::Derivative{1})
    (; ts, cderivs, Xoffset,) = f
    (; i,) = node
    cs = cderivs[1]  # spline coefficients associated to first derivative
    t = ntuple(j -> ts[i - 2 + j], Val(3))  # = (ts[i - 1], ts[i], ts[i + 1])
    y = (
        (t[3] - t[2]) * cs[i - 1] +
        (t[2] - t[1]) * cs[i]
    ) / (t[3] - t[1])
    deperiodise_spline(y, Xoffset, ts, t[2], Val(1))  # only useful if Xoffset ≠ 0 ("infinite" / non-closed filaments)
end

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

# Insertion is based on the knot insertion algorithm for splines.
function insert_node!(f::ClosedSplineFilament, i::Integer, ζ::Real)
    (; ts, cs, Xs,) = f
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]  # insert knot somewhere in the middle (ζ ∈ [0, 1]; by default ζ = 0.5)
    spline_insert_knot!(cs, ts, i, t)    # note: this calls pad_periodic! on cs and ts
    Xnew = f(t; ileft = i)
    insert!(Xs, i + 1, Xnew)
    Xnew
end

function remove_node!(f::ClosedSplineFilament, i::Integer)
    (; ts, cs, Xs,) = f
    T = ts[end + 1] - ts[begin]
    popat!(cs, i)
    popat!(ts, i)
    popat!(Xs, i)
    # Avoid error if number of nodes decreased below the limit.
    # Note that `check_nodes` uses the length of f.Xs to determine this.
    if check_nodes(Bool, f)
        pad_periodic!(cs)
        pad_periodic!(ts, T)
    end
    nothing
end

# If knots were not removed but only inserted, one can pass `removed = false` to avoid
# reevaluating node positions (not needed for node insertion, since the curve is unchanged
# in that case).
function update_after_changing_nodes!(f::ClosedSplineFilament; removed = true)
    (; Xs, cs, ts,) = f
    @assert length(f) == length(Xs) == length(cs) == length(ts)  # these already have the right size
    resize!(f, length(Xs))   # resize all vectors in the filament
    check_nodes(Bool, f) || return f  # avoids error if the new number of nodes is too low
    pad_periodic!(Xs, f.Xoffset)
    # If we removed nodes, recompute spline coefficients (i.e. re-interpolate from remaining nodes).
    _update_coefficients_only!(f; only_derivatives = !removed)
    f
end
