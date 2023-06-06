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
end

function ClosedSplineFilament(N::Integer, ::Type{T}; offset = zero(Vec3{T})) where {T}
    M = 3  # padding needed for cubic splines
    ts = PaddedVector{M}(Vector{T}(undef, N + 2M))
    Xs = similar(ts, Vec3{T})
    cs = similar(Xs)
    ċs = similar(Xs)
    c̈s = similar(Xs)
    cderivs = (ċs, c̈s)
    Xoffset = convert(Vec3{T}, offset)
    ClosedSplineFilament(ts, Xs, cs, cderivs, Xoffset)
end

function change_offset(f::ClosedSplineFilament{T}, offset::Vec3) where {T}
    Xoffset = convert(Vec3{T}, offset)
    ClosedSplineFilament(f.ts, f.Xs, f.cs, f.cderivs, Xoffset)
end

allvectors(f::ClosedSplineFilament) = (f.ts, f.Xs, f.cs, f.cderivs...)

init(::Type{ClosedFilament{T}}, N::Integer, ::CubicSplineMethod; kws...) where {T} =
    ClosedSplineFilament(N, T; kws...)

function Base.similar(f::ClosedSplineFilament, ::Type{T}, dims::Dims{1}) where {T <: Number}
    N, = dims
    ClosedSplineFilament(N, T; offset = f.Xoffset)
end

discretisation_method(::ClosedSplineFilament) = CubicSplineMethod()
interpolation_method(::ClosedSplineFilament) = CubicSplineMethod()

function _update_coefficients_only!(f::ClosedSplineFilament; only_derivatives = false)
    (; ts, Xs, cs, cderivs, Xoffset,) = f
    if !only_derivatives
        solve_cubic_spline_coefficients!(cs, ts, Xs; buf = cderivs[1], Xoffset,)
    end
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

# Refinement is based on the knot insertion algorithm for splines.
function refine!(f::ClosedSplineFilament, crit::RefinementCriterion)
    (; ts, cs, Xs,) = f
    N = length(cs)  # original number of nodes

    # Determine indices of nodes to modify (refine or remove).
    cache = _nodes_to_refine!(f, crit)
    (; inds, remove,) = cache
    n_modify = length(inds)
    iszero(n_modify) && return (n_modify, n_modify)  # = (n_add = 0, n_rem = 0)

    n_rem = sum(remove)  # note: `remove` is a vector of Bool
    n_add = n_modify - n_rem
    @assert n_add ≥ 0

    # Worst case scenario: we add all knots first, then we remove all knots to be removed.
    if n_add > 0
        sizehint!(cs, N + n_add)
        sizehint!(ts, N + n_add)
    end

    # We iterate in reverse to avoiding the need to shift indices.
    for n ∈ reverse(eachindex(inds))
        i, rem = inds[n], remove[n]
        if rem
            popat!(cs, i)
            popat!(ts, i)
            popat!(Xs, i)
        else
            t = (ts[i] + ts[i + 1]) / 2        # insert knot in the middle of the segment
            spline_insert_knot!(cs, ts, i, t)  # note: this calls pad_periodic! on cs and ts
            insert!(Xs, i + 1, f(t; ileft = i))
        end
    end

    if n_add + n_rem > 0
        @assert length(cs) == N + n_add - n_rem
        resize!(f, length(cs))   # resize all vectors in the filament
        if check_nodes(Bool, f)  # avoids error if the new number of nodes is too low
            pad_periodic!(Xs, f.Xoffset)  # just in case...
            _update_coefficients_only!(f; only_derivatives = true)
        end
    end

    n_add, n_rem
end

# Coefficients should be updated before refinement, but not after (since we use
# knot insertion + updating of derivatives).
update_coefficients_before_refinement(::ClosedSplineFilament) = true
update_coefficients_after_refinement(::ClosedSplineFilament) = false
