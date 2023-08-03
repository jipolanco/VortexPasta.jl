export
    ClosedLocalFilament

@doc raw"""
    ClosedLocalFilament{T, M} <: ClosedFilament{T} <: AbstractFilament{T} <: AbstractVector{Vec3{T}}

Describes a closed curve (a loop) in 3D space using a local discretisation
method such as [`FiniteDiffMethod`](@ref).

The parameter `M` corresponds to the number of neighbouring data points needed
on each side of a given discretisation point to estimate derivatives.

---

    Filaments.init(
        ClosedFilament{T}, N::Integer, discret::LocalDiscretisationMethod;
        offset = zero(Vec3{T}),
    ) -> ClosedLocalFilament{T}

Allocate data for a closed filament with `N` discretisation points.

The element type `T` can be omitted, in which case the default `T = Float64` is used.

See also [`Filaments.init`](@ref).

# Examples

Initialise filament with set of discretisation points:

```jldoctest ClosedLocalFilament
julia> f = Filaments.init(ClosedFilament, 16, FiniteDiffMethod(2, HermiteInterpolation(2)));

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

```jldoctest ClosedLocalFilament
julia> f[4, Derivative(1)]
3-element SVector{3, Float64} with indices SOneTo(3):
  0.975687843883729
 -0.7190396901563083
  0.0

julia> f[4, Derivative(2)]
3-element SVector{3, Float64} with indices SOneTo(3):
  0.037089367352557384
 -0.5360868773346441
  0.0
```

Estimate coordinates and derivatives in-between discretisation points:

```jldoctest ClosedLocalFilament
julia> f(4, 0.32)
3-element SVector{3, Float64} with indices SOneTo(3):
 -0.15995621299009463
 -1.1254976317779821
  0.0

julia> Ẋ, Ẍ = f(4, 0.32, Derivative(1)), f(4, 0.32, Derivative(2))
([0.8866970267571707, -0.9868145656366478, 0.0], [-0.6552471551692449, -0.5406810630674177, 0.0])

julia> X′, X″ = normalise_derivatives(Ẋ, Ẍ)
([0.6683664614205478, -0.7438321539488431, 0.0], [-0.358708958397356, -0.32231604392350593, 0.0])

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
struct ClosedLocalFilament{
        T <: AbstractFloat,
        M,  # padding
        Discretisation <: LocalDiscretisationMethod,
        Knots <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
    } <: ClosedFilament{T}

    discretisation :: Discretisation  # derivative estimation method

    # Parametrisation knots: t_i = ∑_{j = 1}^{i - 1} ℓ_i
    ts      :: Knots

    # Discretisation points X_i.
    Xs      :: Points

    # Derivatives (∂ₜX)_i and (∂ₜₜX)_i at nodes with respect to the parametrisation `t`.
    Xderivs :: NTuple{2, Points}

    Xoffset :: Vec3{T}
end

function ClosedLocalFilament(
        N::Integer, discretisation::LocalDiscretisationMethod, ::Type{T};
        offset = zero(Vec3{T}),
    ) where {T}
    M = npad(discretisation)
    ts = PaddedVector{M}(Vector{T}(undef, N + 2M))
    Xs = similar(ts, Vec3{T})
    Xderivs = (similar(Xs), similar(Xs))
    Xoffset = convert(Vec3{T}, offset)
    ClosedLocalFilament(discretisation, ts, Xs, Xderivs, Xoffset)
end

function change_offset(f::ClosedLocalFilament{T}, offset::Vec3) where {T}
    Xoffset = convert(Vec3{T}, offset)
    ClosedLocalFilament(f.discretisation, f.ts, f.Xs, f.Xderivs, Xoffset)
end

allvectors(f::ClosedLocalFilament) = (f.ts, f.Xs, f.Xderivs...)

function init(
        ::Type{ClosedFilament{T}}, N::Integer, disc::LocalDiscretisationMethod;
        kws...,
    ) where {T}
    ClosedLocalFilament(N, disc, T; kws...)
end

function Base.similar(f::ClosedLocalFilament, ::Type{T}, dims::Dims{1}) where {T <: Number}
    N, = dims
    ClosedLocalFilament(N, discretisation_method(f), T; offset = f.Xoffset)
end

discretisation_method(f::ClosedLocalFilament) = f.discretisation
interpolation_method(f::ClosedLocalFilament) = interpolation_method(discretisation_method(f))

function _update_coefficients_only!(f::ClosedLocalFilament)
    (; ts, Xs, Xderivs,) = f
    M = npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts
    method = discretisation_method(f)
    @inbounds for i ∈ eachindex(Xs)
        ℓs_i = ntuple(j -> ts[i - M + j] - ts[i - M - 1 + j], Val(2M))  # = ℓs[(i - M):(i + M - 1)]
        Xs_i = ntuple(j -> Xs[i - M - 1 + j], Val(2M + 1))  # = Xs[(i - M):(i + M)]
        coefs_dot = coefs_first_derivative(method, ℓs_i)
        coefs_ddot = coefs_second_derivative(method, ℓs_i)
        Xderivs[1][i] = sum(splat(*), zip(coefs_dot, Xs_i))  # = ∑ⱼ c[j] * x⃗[j]
        Xderivs[2][i] = sum(splat(*), zip(coefs_ddot, Xs_i))
    end
    # These paddings are needed for Hermite interpolations and stuff like that.
    # (In principle we just need M = 1 for two-point Hermite interpolations.)
    map(pad_periodic!, Xderivs)
    f
end

function (f::ClosedLocalFilament)(node::AtNode, ::Derivative{n} = Derivative(0)) where {n}
    (; Xs, Xderivs,) = f
    coefs = (Xs, Xderivs...)
    coefs[n + 1][node.i]
end

function (f::ClosedLocalFilament)(i::Int, ζ::Number, deriv::Derivative = Derivative(0))
    m = interpolation_method(f)
    _interpolate(m, f, i, ζ, deriv)  # ζ should be in [0, 1]
end

function (f::ClosedLocalFilament)(
        t_in::Number, deriv::Derivative = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    )
    (; ts,) = f
    i, t = _find_knot_segment(ileft, knotlims(f), ts, t_in)
    ζ = (t - ts[i]) / (ts[i + 1] - ts[i])
    y = f(i, ζ, deriv)
    _deperiodise_finitediff(deriv, y, f, t, t_in)
end

function _deperiodise_finitediff(::Derivative{0}, y, f::ClosedLocalFilament, t, t_in)
    (; Xoffset,) = f
    Xoffset === zero(Xoffset) && return y
    t == t_in && return y
    ta, tb = knotlims(f)
    T = tb - ta
    dt = t_in - t  # expected to be a multiple of T
    y + dt / T * Xoffset
end

_deperiodise_finitediff(::Derivative, y, args...) = y  # derivatives (n ≥ 1): shift not needed

function _interpolate(
        method::HermiteInterpolation{M}, f::ClosedLocalFilament,
        i::Int, t::Number, deriv::Derivative{N},
    ) where {M, N}
    (; ts, Xs, Xderivs,) = f
    checkbounds(f, i)
    @assert npad(Xs) ≥ 1
    @assert all(X -> npad(X) ≥ 1, Xderivs)
    @inbounds ℓ_i = ts[i + 1] - ts[i]
    @assert M ≤ 2
    α = 1 / (ℓ_i^N)  # normalise returned derivative by ℓ^N
    values_full = (Xs, Xderivs...)
    ℓ_norm = Ref(one(ℓ_i))
    values_i = ntuple(Val(M + 1)) do m
        @inline
        X = values_full[m]
        @inbounds data = (X[i], X[i + 1]) .* ℓ_norm  # normalise m-th derivative by ℓ^m
        ℓ_norm[] *= ℓ_i
        data
    end
    α * interpolate(method, deriv, t, values_i...)
end

# Coefficients should be updated before refinement to make sure we have the right curvatures
# and knots consistent with segment lengths. After refinement it's not needed, since we call
# update_coefficients! in refine! (see above).
update_coefficients_before_refinement(::ClosedLocalFilament) = true
update_coefficients_after_refinement(::ClosedLocalFilament) = false
