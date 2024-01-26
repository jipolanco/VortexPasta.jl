export AbstractParametrisation,
       ChordalParametrisation,
       CentripetalParametrisation,
       FourierParametrisation

@doc raw"""
    AbstractParametrisation

Abstract type representing a curve parametrisation.

For a parametric curve ``\bm{s}(t)``, the choice of the discrete knots ``t_i`` at given
discrete nodes ``\bm{s}_i`` is not unique.
The parametrisation determines the way the knots ``t_i`` are chosen.

Some predefined parametrisations include:

- [`ChordalParametrisation`](@ref),
- [`CentripetalParametrisation`](@ref),
- [`FourierParametrisation`](@ref).

Moreover, [`CustomParametrisation`](@ref) allows to define custom parametrisation functions.
"""
abstract type AbstractParametrisation <: Function end

@doc raw"""
    ChordalParametrisation <: AbstractParametrisation
    ChordalParametrisation()

Represents the *chordal* parametrisation.

In this case, the increment between two knots is roughly equal to the length of the segment
between the knots:

```math
t_{i + 1} = t_{i} + |\bm{s}_{i + 1} - \bm{s}_i|
```

This is known as the chordal parametrisation in the context of [Catmull–Rom
splines](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline).

This is the recommended (and default) parametrisation for general cases.
"""
struct ChordalParametrisation <: AbstractParametrisation end
(::ChordalParametrisation)(Xs, i) = @inbounds @fastmath norm(Xs[i + 1] - Xs[i])

@doc raw"""
    CentripetalParametrisation <: AbstractParametrisation
    CentripetalParametrisation()

Represents the *centripetal* parametrisation.

In this case, the increment between two knots is given by:

```math
t_{i + 1} = t_{i} + |\bm{s}_{i + 1} - \bm{s}_i|^{1/2}
```

This is known as the centripetal parametrisation in the context of [Catmull–Rom
splines](https://en.wikipedia.org/wiki/Centripetal_Catmull%E2%80%93Rom_spline).
"""
struct CentripetalParametrisation <: AbstractParametrisation end
(::CentripetalParametrisation)(Xs, i) = @inbounds @fastmath sqrt(norm(Xs[i + 1] - Xs[i]))

# Constant parametrisation with 2π period.
# This is *required* by FourierMethod! (Assumed in particular when estimating derivatives.)
@doc raw"""
    FourierParametrisation <: AbstractParametrisation
    FourierParametrisation()
    
Curve parametrisation adapted for Fourier series representation of curves.

In this case, the increment between two knots is constant and given by:

```math
t_{i + 1} = t_{i} + \frac{2π}{N}
```

where ``N`` is the number of nodes.
"""
struct FourierParametrisation <: AbstractParametrisation end
(::FourierParametrisation)(Xs, i) = 2π / length(Xs)

@doc raw"""
    CustomParametrisation{F <: Function} <: AbstractParametrisation
    CustomParametrisation(f::Function)

Allows to define a custom curve parametrisation.

Here the function `f` should have the signature `f(Xs, i) → Δt` where `Xs` is the vector of
discretisation nodes. It must return the knot increment `Δt = t_{i + 1} - t_{i}`.

For instance, the [`ChordalParametrisation`](@ref) corresponds to defining

    f(Xs, i) = norm(Xs[i + 1] - Xs[i])

"""
struct CustomParametrisation{F <: Function} <: AbstractParametrisation
    f::F
end
(p::CustomParametrisation)(Xs, i) = p.f(Xs, i)

to_parametrisation(p::AbstractParametrisation) = p
to_parametrisation(f::F) where {F <: Function} = CustomParametrisation(f)
