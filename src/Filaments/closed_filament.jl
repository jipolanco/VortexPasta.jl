export ClosedFilament

@doc raw"""
    ClosedFilament{T, D <: DiscretisationMethod} <: AbstractFilament{T}

Describes a closed curve (a loop) in 3D space.

It can also be used to represent infinite but unclosed curves described by a periodic function
(such as an infinite straight line or a sinusoidal curve).

`ClosedFilament`s should be generally constructed using [`Filaments.init`](@ref).

# Extended help

## Examples

The following examples use the [`CubicSplineMethod`](@ref) for representing filament curves,
but other methods are also available.

Initialise filament from a set of discretisation points:

```jldoctest ClosedFilament; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
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

```jldoctest ClosedFilament; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
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

```jldoctest ClosedFilament; filter = r"(\d*)\.(\d{13})\d+" => s"\1.\2***"
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

## Curve parametrisation

The parametrisation knots ``t_i`` are directly obtained from the interpolation point
positions.
A standard choice, which is used here by default, is for the knot increments to
approximate the arc length between two interpolation points:

```math
ℓ_{i} ≡ t_{i + 1} - t_{i} = |\bm{s}_{i + 1} - \bm{s}_i|,
```

which is a zero-th order approximation (and a lower bound) for the actual
arc length between points ``\bm{s}_i`` and ``\bm{s}_{i + 1}``.

"""
struct ClosedFilament{
        T,
        Method <: DiscretisationMethod,
        Parametrisation <: AbstractParametrisation,
        M,  # padding (for dealing with periodicity)
        Knots <: PaddedVector{M, T},
        Points <: PaddedVector{M, Vec3{T}},
        Coefs <: DiscretisationCoefs{Vec3{T}, Method},
    } <: AbstractFilament{T, Method, Parametrisation}

    parametrisation :: Parametrisation  # parametrisation function

    # Parametrisation knots tᵢ.
    ts :: Knots

    # Discretisation points s⃗ᵢ.
    Xs :: Points

    # Discretisation and interpolation coefficients.
    coefs :: Coefs

    # End-to-end offset.
    Xoffset :: Vec3{T}
end

function ClosedFilament(
        parametrisation::F, Xs::PaddedVector{M, Vec3{T}}, method::DiscretisationMethod;
        offset = zero(Vec3{T}), nderivs::Val{Nderivs} = Val(2),
    ) where {F <: Function, M, T, Nderivs}
    @assert M == npad(method)
    ts = similar(Xs, T)
    coefs = init_coefficients(method, Xs, nderivs)
    C = continuity(method)
    # Always allow Nderivs ≤ 2
    max(C, 2) ≥ Nderivs || throw(ArgumentError(lazy"$method only allows up to $C derivatives"))
    Xoffset = convert(Vec3{T}, offset)
    ClosedFilament(to_parametrisation(parametrisation), ts, Xs, coefs, Xoffset)
end

ClosedFilament(Xs::PaddedVector, method::DiscretisationMethod; kws...) =
    ClosedFilament(default_parametrisation(method), Xs, method; kws...)

discretisation_method(::Type{<:ClosedFilament{T, D}}) where {T, D} = D()
discretisation_method(f::ClosedFilament) = discretisation_method(typeof(f))

interpolation_method(f::ClosedFilament) = interpolation_method(discretisation_method(f))

init(::Type{ClosedFilament}, N::Integer, args...; kws...) =
    init(ClosedFilament{Float64}, N, args...; kws...)

init(func::Func, ::Type{ClosedFilament}, args...; kws...) where {Func} =
    init(func, ClosedFilament{Float64}, args...; kws...)

function init(
        ::Type{ClosedFilament{T}}, N::Integer, method::DiscretisationMethod;
        parametrisation::F = default_parametrisation(method), kws...,
    ) where {T, F}
    M = npad(method)
    Xs = PaddedVector{M}(Vector{Vec3{T}}(undef, N + 2M))
    ClosedFilament(parametrisation, Xs, method; kws...)
end

function init(
        ::Type{<:ClosedFilament}, Xs::PaddedVector{M, <:Vec3}, method::DiscretisationMethod;
        parametrisation::F = default_parametrisation(method), kws...,
    ) where {M, F}
    @assert M == npad(method)
    f = ClosedFilament(parametrisation, Xs, method; kws...)
    update_coefficients!(f)
end

function init(
        S::Func, ::Type{ClosedFilament{T}},
        τs::AbstractVector, method::DiscretisationMethod;
        kws...,
    ) where {Func <: Function, T}
    offset = S(1) .- S(0)
    f = init(ClosedFilament{T}, length(τs), method; offset, kws...)
    @assert eachindex(f) == eachindex(τs)
    for (i, τ) ∈ pairs(τs)
        f[i] = S(τ)
    end
    update_coefficients!(f)
    f
end

function init(
        S::Func, ::Type{ClosedFilament{T}}, N::Integer, args...;
        kws...,
    ) where {Func, T}
    τs = range(0, 1; length = 2N + 1)[2:2:2N]
    init(S, ClosedFilament{T}, τs, args...; kws...)
end

function knotlims(f::ClosedFilament)
    ts = knots(f) :: PaddedVector
    ts[begin], ts[end + 1]  # we can index at end+1 thanks to padding
end

"""
    end_to_end_offset(f::ClosedFilament{T}) -> Vec3{T}

Return the end-to-end offset `Δ⃗ = f[end + 1] - f[begin]` of a "closed" filament.

For actually closed filaments, the end-to-end offset is zero. However, `ClosedFilament` also
supports the case of infinite (but unclosed) filaments, which infinitely extend along one or
more Cartesian directions. The restriction imposed by `ClosedFilament` is that infinite
filaments repeat themselves periodically, such that `f[i + m * N] == f[i] + m * Δ⃗` where `N`
is the `length` of the filament (i.e. the number of degrees of freedom, or the total number
of *independent* filament nodes).
"""
end_to_end_offset(f::ClosedFilament) = f.Xoffset

"""
    change_offset(f::ClosedFilament, offset::Vec3{<:Real}) -> f′

Change spatial offset between a filament start and endpoints.

See [`init`](@ref) for more details on optional spatial offsets.

This function is allocation-free. It returns a new filament which shares the same
arrays as `f`, and only differs in the offset. Modifying nodes of the returned
filament also modifies nodes of `f`.
"""
function change_offset(f::ClosedFilament{T}, offset::Vec3) where {T}
    Xoffset = convert(Vec3{T}, offset)
    ClosedFilament(f.parametrisation, f.ts, f.Xs, f.coefs, Xoffset)
end

allvectors(f::ClosedFilament) = (f.ts, f.Xs, allvectors(f.coefs)...)

function Base.similar(
        f::ClosedFilament, ::Type{T}, dims::Dims{1};
        offset = f.Xoffset,
        method = discretisation_method(f),
        nderivs::Val = Val(nderivatives(f)),
    ) where {T <: Number}
    Xs = similar(nodes(f), Vec3{T}, dims)
    ClosedFilament(f.parametrisation, Xs, method; offset, nderivs)
end

function update_coefficients!(f::ClosedFilament; knots = nothing)
    (; parametrisation, ts, Xs,) = f
    Xoffset = end_to_end_offset(f)
    M = npad(Xs)
    check_nodes(f)

    # 1. Periodically pad Xs.
    pad_periodic!(Xs, Xoffset)

    # 2. Compute parametrisation knots `ts`.
    @assert M == npad(ts)
    @assert M ≥ 1  # minimum padding required for computation of ts
    _update_knots_periodic!(parametrisation, ts, Xs, knots)

    # 3. Estimate coefficients needed for derivatives and interpolations.
    _update_coefficients_only!(f)

    f
end

function _update_coefficients_only!(f::ClosedFilament; kws...)
    # Dispatch to the implementation associated to the chosen discretisation method.
    method = discretisation_method(f)
    _update_coefficients_only!(method, f; kws...)
end

## Evaluation of values and derivatives on discretisation points (nodes)
(f::ClosedFilament)(node::AtNode, ::Derivative{0} = Derivative(0)) = f[node.i]

# We let the different discretisation methods provide (possibly optimised) implementations
# for derivatives at nodes:
(f::ClosedFilament)(node::AtNode, d::Derivative) =
    _derivative_at_node(d, discretisation_method(f), f::ClosedFilament, node::AtNode)

## Interpolation of values and derivatives in-between nodes

# 1. If we know the interpolation segment (here ζ ∈ [0, 1])
function (f::ClosedFilament{T})(i::Int, ζ_in::Number, d::Derivative = Derivative(0)) where {T}
    ζ = _convert_float(T, ζ_in)
    _interpolate(discretisation_method(f), f, i, ζ, d)
end

# 2. If we don't know the interpolation segment (here `t` is the curve parameter).
function (f::ClosedFilament{T})(
        t_in::Number, d::Derivative = Derivative(0);
        ileft::Union{Nothing, Int} = nothing,
    ) where {T}
    t = _convert_float(T, t_in)
    _interpolate(discretisation_method(f), f, t, d; ileft)
end

# Convert value to type T. We don't perform conversion if the value is not a float.
# This is to workaround issue when working with ForwardDiff.Dual (which is never an
# AbstractFloat).
_convert_float(::Type{T}, x::AbstractFloat) where {T} = convert(T, x)
_convert_float(::Type{T}, x) where {T} = x

## Knot insertion and removal.
#  We allow the different discretisation methods to implement their own versions (for
#  example, using standard knot insertion algorithms for splines).

insert_node!(f::ClosedFilament, i::Integer, ζ::Real) =
    _insert_node!(discretisation_method(f), f, i, ζ)

remove_node!(f::ClosedFilament, i::Integer) =
    _remove_node!(discretisation_method(f), f, i)

update_after_changing_nodes!(f::ClosedFilament; removed = true) =
    _update_after_changing_nodes!(discretisation_method(f), f; removed)

## Default implementation

function _insert_node!(::DiscretisationMethod, f::ClosedFilament, i::Integer, ζ::Real)
    (; Xs, ts,) = f
    (; cs, cderivs,) = f.coefs
    nderivs = length(cderivs)

    # New point to be added.
    # Note: this uses derivatives at nodes (Hermite interpolations), so make sure they are up to date!
    s⃗ = f(i, ζ)

    # We also insert new derivatives in case we need to perform Hermite interpolations later
    # (for instance if we insert other nodes).
    # Derivatives will be recomputed later when calling `update_after_changing_nodes!`.
    derivs = _eval_derivatives_up_to(Val(nderivs), f, i, ζ)

    # New knot location
    t = (1 - ζ) * ts[i] + ζ * ts[i + 1]

    insert!(ts, i + 1, t)
    insert!(Xs, i + 1, s⃗)
    insert!(cs, i + 1, s⃗)
    for (coefs, val) ∈ zip(cderivs, derivs)
        insert!(coefs, i + 1, val)
    end

    s⃗
end

@inline function _eval_derivatives_up_to(::Val{n}, f, i, ζ) where {n}
    prev = _eval_derivatives_up_to(Val(n - 1), f, i, ζ)
    val = f(i, ζ, Derivative(n))
    (prev..., val)
end

@inline _eval_derivatives_up_to(::Val{0}, args...) = ()

function _remove_node!(::DiscretisationMethod, f::ClosedFilament, i::Integer)
    (; ts, Xs,) = f
    # No need to modify interpolation coefficients, since they will be updated later in
    # `update_after_changing_nodes!`. This assumes that we won't insert nodes after removing
    # nodes (i.e. insertions are done before removals).
    popat!(ts, i)
    popat!(Xs, i)
end

# The `removed` argument is just there for compatibility with splines.
function _update_after_changing_nodes!(::DiscretisationMethod, f::ClosedFilament; removed = true)
    (; Xs,) = f
    resize!(f, length(Xs))   # resize all vectors in the filament
    if check_nodes(Bool, f)  # avoids error if the new number of nodes is too low
        pad_periodic!(Xs, f.Xoffset)
        update_coefficients!(f)
    end
    f
end
