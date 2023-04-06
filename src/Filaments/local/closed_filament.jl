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

# Interpolating in-between discretisation points

One can estimate filament coordinates and derivatives in-between discretisation
points using the chosen interpolation method.
The signature for doing this is as follows:

    f = ClosedLocalFilament(...)
    f(i::Int, t::Number, [derivative = Derivative(0)])

This returns an interpolation between nodes ``\bm{X}_i`` and ``\bm{X}_{i + 1}``.
Here, the ``t`` parameter should be in ``[0, 1]``.
To obtain `N`-th order derivatives, pass `Derivative(N)` as the last argument.

Note that, if using Hermite interpolations (which is the default), one must
first estimate filament derivatives using [`update_coefficients!`](@ref).

# Examples

```jldoctest
julia> fil = Filaments.init(ClosedFilament, 16, FiniteDiffMethod(2); interpolation_method = HermiteInterpolation(2));

julia> θs = range(-1, 1; length = 17)[1:16]
-1.0:0.125:0.875

julia> @. fil = Vec3(cospi(θs), sinpi(θs), 0)
16-element ClosedLocalFilament{Float64, FiniteDiffMethod{2}}:
 [-1.0, -0.0, 0.0]
 [-0.9238795325112867, -0.3826834323650898, 0.0]
 [-0.7071067811865476, -0.7071067811865476, 0.0]
 [-0.3826834323650898, -0.9238795325112867, 0.0]
 [0.0, -1.0, 0.0]
 [0.3826834323650898, -0.9238795325112867, 0.0]
 [0.7071067811865476, -0.7071067811865476, 0.0]
 [0.9238795325112867, -0.3826834323650898, 0.0]
 [1.0, 0.0, 0.0]
 [0.9238795325112867, 0.3826834323650898, 0.0]
 [0.7071067811865476, 0.7071067811865476, 0.0]
 [0.3826834323650898, 0.9238795325112867, 0.0]
 [0.0, 1.0, 0.0]
 [-0.3826834323650898, 0.9238795325112867, 0.0]
 [-0.7071067811865476, 0.7071067811865476, 0.0]
 [-0.9238795325112867, 0.3826834323650898, 0.0]

julia> fil[4]
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.3826834323650898
 -0.9238795325112867
  0.0

julia> fil[5] = (fil[4] + 2 * fil[6]) ./ 2
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
  0.1913417161825449
 -1.38581929876693
  0.0

julia> update_coefficients!(fil);

julia> derivative(fil, 1)
16-element PaddedVector{2, StaticArraysCore.SVector{3, Float64}, Vector{StaticArraysCore.SVector{3, Float64}}}:
 [-1.1102230246251565e-16, -1.0056712250866777, 0.0]
 [0.38485371624697473, -0.9291190612931332, 0.0]
 [0.7182731830606435, -0.6888706857208597, 0.0]
 [0.9756878438837288, -0.7190396901563081, 0.0]
 [0.5108238194008445, 0.34240719514000695, 0.0]
 [0.6420265791563875, 0.7226786074475605, 0.0]
 [0.7305563783649631, 0.6383040901696273, 0.0]
 [0.3848537162469746, 0.9291190612931333, 0.0]
 [1.1102230246251565e-16, 1.0056712250866777, 0.0]
 [-0.38485371624697473, 0.9291190612931332, 0.0]
 [-0.7111169429029729, 0.7111169429029729, 0.0]
 [-0.9291190612931333, 0.3848537162469746, 0.0]
 [-1.0056712250866777, 1.1102230246251565e-16, 0.0]
 [-0.9291190612931332, -0.38485371624697473, 0.0]
 [-0.7111169429029729, -0.7111169429029729, 0.0]
 [-0.3848537162469746, -0.9291190612931333, 0.0]

julia> fil(4, 0.32)
3-element StaticArraysCore.SVector{3, Float64} with indices SOneTo(3):
 -0.1599562129900946
 -1.1254976317779821
  0.0

julia> Ẋ, Ẍ = fil(4, 0.32, Derivative(1)), fil(4, 0.32, Derivative(2))
([0.8866970267571704, -0.9868145656366478, 0.0], [-0.6552471551692444, -0.5406810630674178, 0.0])

julia> X′, X″ = normalise_derivatives(Ẋ, Ẍ)
([0.6683664614205477, -0.7438321539488432, 0.0], [-0.358708958397356, -0.3223160439235058, 0.0])

```

---

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
        Distances <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
    } <: ClosedFilament{T}

    discretisation :: Discretisation  # derivative estimation method
    interpolation  :: Interpolation   # interpolation method

    # Node distances ℓ_i = |X_{i + 1} - X_i|.
    ℓs      :: Distances

    # Discretisation points X_i.
    Xs      :: Points

    # First derivatives (∂ₜX)_i at nodes with respect to the parametrisation `t`.
    Xs_dot  :: Points

    # Second derivatives (∂ₜₜX)_i at nodes with respect to the parametrisation `t`.
    Xs_ddot :: Points

    function ClosedLocalFilament(
            N::Integer, discretisation::LocalDiscretisationMethod,
            interpolation::LocalInterpolationMethod, ::Type{T},
        ) where {T}
        M = npad(discretisation)
        ℓs = PaddedVector{M}(Vector{T}(undef, N + 2M))
        Xs = similar(ℓs, Vec3{T})
        Xs_dot = similar(Xs)
        Xs_ddot = similar(Xs)
        new{T, M, typeof(discretisation), typeof(interpolation), typeof(ℓs), typeof(Xs)}(
            discretisation, interpolation, ℓs, Xs, Xs_dot, Xs_ddot,
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
derivatives(f::ClosedLocalFilament) = (f.Xs_dot, f.Xs_ddot)

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
    (; ℓs, Xs, Xs_dot, Xs_ddot,) = f

    # 1. Periodically pad Xs.
    pad_periodic!(Xs)

    # 2. Compute node distances ℓs.
    M = npad(Xs)
    @assert M == npad(ℓs)
    @assert M ≥ 1  # minimum padding required for computation of ℓs
    @assert eachindex(Xs) == eachindex(ℓs)
    @inbounds for i ∈ eachindex(ℓs)
        ℓs[i] = norm(Xs[i + 1] - Xs[i])
    end
    pad_periodic!(ℓs)

    # 3. Estimate derivatives at nodes.
    method = discretisation_method(f)
    @inbounds for i ∈ eachindex(Xs_dot)
        ℓs_i = ntuple(j -> ℓs[i - M - 1 + j], Val(2M))      # = ℓs[(i - M):(i + M - 1)]
        Xs_i = ntuple(j -> Xs[i - M - 1 + j], Val(2M + 1))  # = Xs[(i - M):(i + M)]
        coefs_dot = coefs_first_derivative(method, ℓs_i)
        coefs_ddot = coefs_second_derivative(method, ℓs_i)
        Xs_dot[i] = sum(splat(*), zip(coefs_dot, Xs_i))  # = ∑ⱼ c[j] * x⃗[j]
        Xs_ddot[i] = sum(splat(*), zip(coefs_ddot, Xs_i))
    end

    # These paddings are needed for Hermite interpolations and stuff like that.
    # (In principle we just need M = 1 for two-point Hermite interpolations.)
    pad_periodic!(Xs_dot)
    pad_periodic!(Xs_ddot)

    f
end

function (f::ClosedLocalFilament)(i::Int, t::Number, deriv::Derivative = Derivative(0))
    m = interpolation_method(f)
    _interpolate(m, f, i, t, deriv)
end

function _interpolate(
        method::HermiteInterpolation{M}, f::ClosedLocalFilament,
        i::Int, t::Number, deriv::Derivative{N},
    ) where {M, N}
    (; ℓs, Xs, Xs_dot, Xs_ddot,) = f
    checkbounds(eachindex(f), i)
    @assert npad(Xs) ≥ 1
    @assert npad(Xs_dot) ≥ 1
    @assert npad(Xs_ddot) ≥ 1
    @inbounds ℓ_i = ℓs[i]
    @assert M ≤ 2
    α = 1 / (ℓ_i^N)  # normalise returned derivative by ℓ^N
    values_full = (Xs, Xs_dot, Xs_ddot)
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
