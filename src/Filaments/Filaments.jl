"""
    Filaments

Module for dealing with the discretisation of curves in 3D space.
"""
module Filaments

export
    ClosedFilament,
    Vec3,
    nodes,
    compute_derivatives!

using Base: @propagate_inbounds
using LinearAlgebra: norm, ⋅, ×
using StaticArrays

"""
    Vec3{T}

Three-element static vector, alias for `SVector{3, T}`.

Used to describe vectors and coordinates in 3D space.
"""
const Vec3{T} = SVector{3, T}

include("derivative_estimation.jl")
include("padded_vector.jl")

@doc raw"""
    ClosedFilament{T, M} <: AbstractVector{Vec3{T}}

Describes a closed curve (a loop) in 3D space.

The curve coordinates are parametrised as ``\bm{X}(t)`` with ``t ∈ ℝ``.

The curve is discretised by a set of *nodes* (or discretisation points)
``\bm{X}(t_i) = \bm{X}_i`` for ``i ∈ \{1, 2, …, N\}``.

A `ClosedFilament` is treated as an `AbstractVector` of length `N`, in which
each element is a discretisation point `\bm{X}_i`. Therefore, one can use the
usual indexing notation to retrieve and to modify discretisation points. See
further below for examples.

---

    ClosedFilament(N::Integer, method::DerivativeEstimationMethod, [T = Float64])

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
Note that this choice is known as *chordal* parametrisation in the context of
[Catmull--Rom splines](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline).
"""
struct ClosedFilament{
        T <: AbstractFloat,
        M,
        DEM <: DerivativeEstimationMethod{M},
        Distances <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
    } <: AbstractVector{Vec3{T}}

    # Derivative estimation method.
    method_derivatives :: DEM

    # Node distances ℓ_i = |X_{i + 1} - X_i|.
    ℓs      :: Distances

    # Discretisation points X_i.
    Xs      :: Points

    # First derivatives (∂ₜX)_i at nodes with respect to the parametrisation `t`.
    Xs_dot  :: Points

    # Second derivatives (∂ₜₜX)_i at nodes with respect to the parametrisation `t`.
    Xs_ddot :: Points

    function ClosedFilament(
            N::Integer, method::DerivativeEstimationMethod, ::Type{T},
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

Returns the discretisation points ``\\bm{X}_i`` of the filament.
"""
nodes(f::ClosedFilament) = f.Xs

derivative_estimation_method(f::ClosedFilament) = f.method_derivatives

Base.eltype(::Type{<:ClosedFilament{T}}) where {T} = T
Base.eltype(f::ClosedFilament) = eltype(typeof(f))
Base.size(f::ClosedFilament) = size(nodes(f))

function Base.showarg(io::IO, f::ClosedFilament, toplevel)
    toplevel || print(io, "::")
    T = eltype(f)
    DEM = typeof(derivative_estimation_method(f))
    print(io, nameof(typeof(f)), '{', T, ',', ' ', DEM, '}')
end

@propagate_inbounds Base.getindex(f::ClosedFilament, i::Int) = nodes(f)[i]
@propagate_inbounds Base.setindex!(f::ClosedFilament, v, i::Int) = nodes(f)[i] = v

function compute_derivatives!(f::ClosedFilament)
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
    method = derivative_estimation_method(f)
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

    nothing
end

end
