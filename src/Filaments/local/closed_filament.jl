export
    ClosedLocalFilament

@doc raw"""
    ClosedLocalFilament{T, M} <: ClosedFilament{T} <: AbstractFilament{T} <: AbstractVector{Vec3{T}}

Describes a closed curve (a loop) in 3D space using a local discretisation
method such as [`FiniteDiffMethod`](@ref).

The parameter `M` corresponds to the number of neighbouring data points needed
on each side of a given discretisation point to estimate derivatives.

---

    Filaments.init(ClosedFilament{T}, N::Integer, discret::LocalDiscretisationMethod, [interp = HermiteInterpolation(2)]) -> ClosedLocalFilament{T}

Allocate data for a closed filament with `N` discretisation points.

The element type `T` can be omitted, in which case the default `T = Float64` is used.

See also [`Filaments.init`](@ref).

# Optional arguments

- `interpolation`: method used to evaluate curve coordinates or derivatives
  in-between discretisation points. Must be a type of [`LocalInterpolationMethod`].
  By default, quintic Hermite interpolation (`HermiteInterpolation(2)`) is used.
  Note that Hermite interpolations use the derivatives estimated via the chosen
  discretisation method.

# Examples

Initialise filament with set of discretisation points:

```jldoctest ClosedLocalFilament
julia> f = Filaments.init(ClosedFilament, 16, FiniteDiffMethod(2), HermiteInterpolation(2));

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

```jldoctest ClosedLocalFilament
julia> f(AtNode(4), Derivative(1))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  0.975687843883729
 -0.7190396901563083
  0.0

julia> f(AtNode(4), Derivative(2))
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  0.037089367352557384
 -0.5360868773346441
  0.0
```

Estimate coordinates and derivatives in-between discretisation points:

```jldoctest ClosedLocalFilament
julia> f(4, 0.32)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
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
approximate the arclength between two interpolation points:

```math
ℓ_{i} ≡ t_{i + 1} - t_{i} = |\bm{X}_{i + 1} - \bm{X}_i|,
```

which is a zero-th order approximation (and a lower bound) for the actual
arclength between points ``\bm{X}_i`` and ``\bm{X}_{i + 1}``.
"""
struct ClosedLocalFilament{
        T <: AbstractFloat,
        M,  # padding
        Discretisation <: LocalDiscretisationMethod,
        Interpolation <: LocalInterpolationMethod,
        Knots <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
    } <: ClosedFilament{T}

    discretisation :: Discretisation  # derivative estimation method
    interpolation  :: Interpolation   # interpolation method

    # Parametrisation knots: t_i = ∑_{j = 1}^{i - 1} ℓ_i
    ts      :: Knots

    # Discretisation points X_i.
    Xs      :: Points

    # Derivatives (∂ₜX)_i and (∂ₜₜX)_i at nodes with respect to the parametrisation `t`.
    Xderivs :: NTuple{2, Points}

    function ClosedLocalFilament(
            N::Integer, discretisation::LocalDiscretisationMethod,
            interpolation::LocalInterpolationMethod, ::Type{T},
        ) where {T}
        M = npad(discretisation)
        ts = PaddedVector{M}(Vector{T}(undef, N + 2M))
        Xs = similar(ts, Vec3{T})
        Xderivs = (similar(Xs), similar(Xs))
        new{T, M, typeof(discretisation), typeof(interpolation), typeof(ts), typeof(Xs)}(
            discretisation, interpolation, ts, Xs, Xderivs,
        )
    end
end

# Use default interpolation method
ClosedLocalFilament(N::Integer, disc::LocalDiscretisationMethod, ::Type{T}) where {T} =
    ClosedLocalFilament(N, disc, HermiteInterpolation(2), T)

init(::Type{ClosedFilament{T}}, N::Integer, disc::LocalDiscretisationMethod, args...) where {T} =
    ClosedLocalFilament(N, disc, args..., T)

"""
    derivatives(f::ClosedLocalFilament) -> (Ẋs, Ẍs)

Return a tuple with the first and second derivatives at the filament nodes.

Each element is a vector with the derivatives estimated at the interpolation points.

This function should be generally called after [`update_coefficients!`](@ref).
"""
derivatives(f::ClosedLocalFilament) = f.Xderivs

"""
    derivative(f::ClosedLocalFilament, i::Int) -> AbstractVector

Return ``i``-th derivative at filament nodes.

Same as `derivatives(f)[i]`. The only allowed values are `i ∈ (1, 2)`.

See [`derivatives`](@ref) for details.
"""
derivative(f::ClosedLocalFilament, i::Int) = derivatives(f)[i]

discretisation_method(f::ClosedLocalFilament) = f.discretisation
interpolation_method(f::ClosedLocalFilament) = f.interpolation

function update_coefficients!(f::ClosedLocalFilament)
    (; ts, Xs, Xderivs,) = f

    # 1. Periodically pad Xs.
    pad_periodic!(Xs)

    # 2. Compute parametrisation knots `ts`.
    M = npad(Xs)
    @assert M == npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts
    _update_knots_periodic!(ts, Xs)

    # 3. Estimate derivatives at nodes.
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
        t::Number, deriv::Derivative = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    )
    (; ts,) = f
    i = if ileft === nothing
        searchsortedlast(ts, t) :: Int
    else
        ileft
    end
    ζ = (t - ts[i]) / (ts[i + 1] - ts[i])
    f(i, ζ, deriv)
end

function _interpolate(
        method::HermiteInterpolation{M}, f::ClosedLocalFilament,
        i::Int, t::Number, deriv::Derivative{N},
    ) where {M, N}
    (; ts, Xs, Xderivs,) = f
    checkbounds(eachindex(f), i)
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
