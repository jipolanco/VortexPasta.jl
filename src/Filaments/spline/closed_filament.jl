export
    ClosedSplineFilament

@doc raw"""
    ClosedSplineFilament{T} <: ClosedFilament{T} <: AbstractFilament{T} <: AbstractVector{Vec3{T}}

Describes a closed curve (a loop) in 3D space using periodic cubic splines.

---

    ClosedSplineFilament(N::Integer, [T = Float64])                    -> ClosedSplineFilament
    Filaments.init(ClosedFilament{T}, N::Integer, CubicSplineMethod()) -> ClosedSplineFilament

Allocate data for a closed spline filament with `N` discretisation points.

See also [`Filaments.init`](@ref).

# Extended help

## Curve parametrisation

The parametrisation knots ``t_i`` are directly obtained from the interpolation point
positions.
A standard choice, which is used here, is for the knot increments to
approximate the arclength between two interpolation points:

```math
ℓ_{i} ≡ t_{i + 1} - t_{i} = |\bm{X}_{i + 1} - \bm{X}_i|,
```

which is a zero-th order approximation (and a lower bound) for the actual
arclength between points ``\bm{X}_i`` and ``\bm{X}_{i + 1}``.
"""
struct ClosedSplineFilament{
        T <: AbstractFloat,
        M,  # padding (fixed to 2?)
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

    function ClosedSplineFilament(N::Integer, ::Type{T}) where {T}
        M = 3  # padding needed for cubic splines
        ts = PaddedVector{M}(Vector{T}(undef, N + 2M))
        Xs = similar(ts, Vec3{T})
        cs = similar(Xs)
        ċs = similar(Xs)
        c̈s = similar(Xs)
        cderivs = (ċs, c̈s)
        new{T, M, typeof(ts), typeof(Xs)}(ts, Xs, cs, cderivs)
    end
end

init(::Type{ClosedFilament{T}}, N::Integer, ::CubicSplineMethod) where {T} =
    ClosedSplineFilament(N, T)

discretisation_method(::ClosedSplineFilament) = CubicSplineMethod()
interpolation_method(::ClosedSplineFilament) = CubicSplineMethod()

# TODO optimise solving for coefficients
# - pass "raw" data only? (instead of PaddedVector's)
function update_coefficients!(f::ClosedSplineFilament)
    (; ts, Xs, cs, cderivs,) = f
    pad_periodic!(Xs)
    _update_knots_periodic!(ts, Xs)
    solve_cubic_spline_coefficients!(cs, ts, Xs; buf = cderivs[1])
    spline_derivative!(cderivs[1], cs, ts, Val(4))
    spline_derivative!(cderivs[2], cderivs[1], ts, Val(3))
    f
end

(f::ClosedSplineFilament)(node::AtNode, ::Derivative{0} = Derivative(0)) = f[node.i]
# TODO can we optimise evaluation of derivatives on nodes? (second derivatives should be particularly easy with cubic splines!)
(f::ClosedSplineFilament)(node::AtNode, der::Derivative) = f(node.i, 0.0, der)

function (f::ClosedSplineFilament)(i::Int, ζ::Number, der::Derivative = Derivative(0))
    (; ts,) = f
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]  # convert ζ ∈ [0, 1] to spline parametrisation
    f(t, der; ileft = i)
end

# Here `t` is in the spline parametrisation.
function (f::ClosedSplineFilament)(
        t::Number, ::Derivative{n} = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    ) where {n}
    (; ts, cs, cderivs,) = f
    i = if ileft === nothing
        searchsortedlast(ts, t) :: Int
    else
        ileft
    end
    coefs = (cs, cderivs...)[n + 1]
    ord = 4 - n
    evaluate_spline(coefs, ts, i, t, Val(ord))
end
