"""
    Filaments

Module for dealing with the discretisation of curves in 3D space.
"""
module Filaments

export
    ClosedFilament,
    Vec3,
    Derivative,
    nodes,
    node_distances,
    estimate_derivatives!,
    normalise_derivatives,
    normalise_derivatives!,
    interpolate,
    derivatives,
    derivative

using Base: @propagate_inbounds
using LinearAlgebra: norm, normalize, ⋅, ×
using StaticArrays
using StructArrays

"""
    Vec3{T}

Three-element static vector, alias for `SVector{3, T}`.

Used to describe vectors and coordinates in 3D space.
"""
const Vec3{T} = SVector{3, T}

"""
    Derivative{N}

Represents the ``N``-th order derivative operator.

Used in particular to interpolate derivatives along filaments.
"""
struct Derivative{N} end
Derivative(N::Int) = Derivative{N}()

@doc raw"""
    AbstractFilament{T} <: AbstractVector{Vec3{T}}

Abstract type representing a curve in 3D space.

The curve coordinates are parametrised as ``\bm{X}(t)`` with ``t ∈ ℝ``.

The curve is discretised by a set of *nodes* (or discretisation points)
``\bm{X}(t_i) = \bm{X}_i`` for ``i ∈ \{1, 2, …, N\}``.

See [`ClosedFilament`](@ref) for a concrete implementation of `AbstractFilament`.

An `AbstractFilament` is treated as an `AbstractVector` of length `N`, in which
each element is a discretisation point `\bm{X}_i`. Therefore, one can use the
usual indexing notation to retrieve and to modify discretisation points. See
[`ClosedFilament`](@ref) for some examples.
"""
abstract type AbstractFilament{T} <: AbstractVector{Vec3{T}} end

include("discretisations.jl")
include("local/padded_vector.jl")
include("local/finitediff.jl")
include("local/interpolation.jl")
include("local/interp_hermite.jl")

@doc raw"""
    ClosedFilament{T, M} <: AbstractFilament{T} <: AbstractVector{Vec3{T}}

Describes a closed curve (a loop) in 3D space.

The parameter `M` corresponds to the number of neighbouring data points needed
on each side of a given discretisation point to estimate derivatives.

---

    ClosedFilament(N::Integer, method::DiscretisationMethod, [T = Float64])

Allocate data for a closed filament with `N` discretisation points.

# Examples

```jldoctest
julia> fil = ClosedFilament(16, FiniteDiff(2), Float64);

julia> θs = range(-1, 1; length = 17)[1:16]
-1.0:0.125:0.875

julia> @. fil = Vec3(cospi(θs), sinpi(θs), 0)
16-element ClosedFilament{Float64, FiniteDiff{2}}:
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

julia> estimate_derivatives!(fil);

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

```

---

# Extended help

## Curve parametrisation

The parametrisation knots ``t_i`` are directly obtained from the interpolation point
positions.
A standard choice is for the knot increments to approximate the arclength
between two interpolation points:

```math
ℓ_{i} ≡ t_{i + 1} - t_{i} = |\bm{X}_{i + 1} - \bm{X}_i|,
```

which is a zero-th order approximation (and a lower bound) for the actual
arclength between points ``\bm{X}_i`` and ``\bm{X}_{i + 1}``.
"""
struct ClosedFilament{
        T <: AbstractFloat,
        M,  # padding
        Disc <: DiscretisationMethod,
        Distances <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
    } <: AbstractFilament{T}

    # Derivative estimation method.
    discretisation :: Disc

    # Node distances ℓ_i = |X_{i + 1} - X_i|.
    ℓs      :: Distances

    # Discretisation points X_i.
    Xs      :: Points

    # First derivatives (∂ₜX)_i at nodes with respect to the parametrisation `t`.
    Xs_dot  :: Points

    # Second derivatives (∂ₜₜX)_i at nodes with respect to the parametrisation `t`.
    Xs_ddot :: Points

    function ClosedFilament(
            N::Integer, method::DiscretisationMethod, ::Type{T},
        ) where {T}
        M = npad(method)
        ℓs = PaddedVector{M}(Vector{T}(undef, N + 2M))
        Xs = PaddedVector{M}(Vector{Vec3{T}}(undef, N + 2M))
        Xs_dot = PaddedVector{M}(Vector{Vec3{T}}(undef, N + 2M))
        Xs_ddot = PaddedVector{M}(Vector{Vec3{T}}(undef, N + 2M))
        new{T, M, typeof(method), typeof(ℓs), typeof(Xs)}(
            method, ℓs, Xs, Xs_dot, Xs_ddot,
        )
    end
end

"""
    nodes(f::ClosedFilament) -> PaddedVector

Return the discretisation points ``\\bm{X}_i`` of the filament.
"""
nodes(f::ClosedFilament) = f.Xs

"""
    node_distances(f::ClosedFilament) -> PaddedVector

Return precomputed node distances ``ℓ_i = |\\bm{X}_{i + 1} - \\bm{X}_i|``.

This function should be generally called after [`estimate_derivatives!`](@ref).
"""
node_distances(f::ClosedFilament) = f.ℓs

"""
    derivatives(f::ClosedFilament) -> (Ẋs, Ẍs)

Return a tuple with the first and second derivatives at the filament nodes.

Each element is a vector with the derivatives estimated at the interpolation points.

This function should be generally called after [`estimate_derivatives!`](@ref).
"""
derivatives(f::ClosedFilament) = (f.Xs_dot, f.Xs_ddot)

"""
    derivative(f::ClosedFilament, i::Int) -> AbstractVector

Return ``i``-th derivative at filament nodes.

Same as `derivatives(f)[i]`. The only allowed values are `i ∈ (1, 2)`.

See [`derivatives`](@ref) for details.
"""
derivative(f::ClosedFilament, i::Int) = derivatives(f)[i]

discretisation_method(f::ClosedFilament) = f.discretisation

Base.eltype(::Type{<:ClosedFilament{T}}) where {T} = T
Base.eltype(f::ClosedFilament) = eltype(typeof(f))
Base.size(f::ClosedFilament) = size(nodes(f))

function Base.showarg(io::IO, f::ClosedFilament, toplevel)
    toplevel || print(io, "::")
    T = eltype(f)
    disc = typeof(discretisation_method(f))
    print(io, nameof(typeof(f)), '{', T, ',', ' ', disc, '}')
end

@propagate_inbounds Base.getindex(f::ClosedFilament, i::Int) = nodes(f)[i]
@propagate_inbounds Base.setindex!(f::ClosedFilament, v, i::Int) = nodes(f)[i] = v

"""
    estimate_derivatives!(f::ClosedFilament) -> (Ẋs, Ẍs)

Estimate first and second derivatives at filament nodes based on the locations
of the discretisation points.

Note that derivatives are with respect to the (arbitrary) parametrisation
``\\bm{X}(t)``, and *not* with respect to the arclength ``ξ = ξ(t)``. In other
words, the returned derivatives do not directly correspond to the unit tangent
and curvature vectors (but they are closely related).

The estimated derivatives are returned by this function as a tuple of vectors.

The derivatives are stored in the `ClosedFilament` object, and can also be
retrieved later by calling [`derivatives`](@ref) or [`derivative`](@ref).
"""
function estimate_derivatives!(f::ClosedFilament)
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

    Xs_dot, Xs_ddot
end

"""
    normalise_derivatives(Ẋ::Vec3, Ẍ::Vec3) -> (X′, X″)
    normalise_derivatives((Ẋ, Ẍ)::NTuple)   -> (X′, X″)

Return derivatives with respect to the arc length ``ξ``, from derivatives with
respect to the parameter ``t``.

The returned derivatives satisfy:

- ``\\bm{X}' ≡ t̂`` is the **unit tangent** vector;

- ``\\bm{X}'' ≡ ρ n̂`` is the **curvature** vector, where ``n̂`` is the normal unit
  vector (with ``t̂ ⋅ n̂ = 0``) and ``ρ = R^{-1}`` is the curvature (and R the
  curvature radius).
"""
function normalise_derivatives(Ẋ::Vec3, Ẍ::Vec3)
    t̂ = normalize(Ẋ)  # unit tangent vector (= X′)
    X″ = (Ẍ - (Ẍ ⋅ t̂) * t̂) ./ sum(abs2, Ẋ)  # curvature vector
    t̂, X″
end

normalise_derivatives(derivs::NTuple{2, Vec3}) = normalise_derivatives(derivs...)

"""
    normalise_derivatives!(fil::AbstractFilament)
    normalise_derivatives!(Ẋ::AbstractVector, Ẍ::AbstractVector)

Normalise derivatives at filament nodes.

Note that filament derivatives are modified, and thus Hermite interpolations
may be incorrect after doing this. If possible, prefer using
[`normalise_derivatives`](@ref), which works on a single filament location at a time.

See [`normalise_derivatives`](@ref) for more details.
"""
function normalise_derivatives!(Ẋ::AbstractVector, Ẍ::AbstractVector)
    derivs = StructArray((Ẋ, Ẍ))
    map!(normalise_derivatives, derivs, derivs)
    (Ẋ, Ẍ)
end

normalise_derivatives!(fil::AbstractFilament) = normalise_derivatives!(derivatives(fil)...)

@doc raw"""
    interpolate(
        method::HermiteInterpolation{M}, f::ClosedFilament, i::Int,
        t::Number, [derivative = Derivative(0)],
    )

Evaluate filament coordinate or derivative at arbitrary point using Hermite interpolation.

The filament is interpolated between nodes ``\bm{X}_i`` and ``\bm{X}_{i + 1}``.

The parameter ``t`` should be normalised to be in ``[0, 1]``.
"""
function interpolate(
        method::HermiteInterpolation{M}, f::ClosedFilament,
        i::Int, t::Number, deriv::Derivative{N} = Derivative(0),
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

end
