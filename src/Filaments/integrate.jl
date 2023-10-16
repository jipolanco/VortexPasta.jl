using Base: @propagate_inbounds

@doc raw"""
    integrate(integrand::Function, f::AbstractFilament, i::Int, quad::AbstractQuadrature)
    integrate(integrand::Function, s::Segment, quad::AbstractQuadrature)

Estimate integral along a filament segment using the chosen quadrature.

Integration is performed in the segment `(f[i], f[i + 1])`.

The function `integrand(ζ)` takes a single argument ``ζ ∈ [0, 1]`` which corresponds
to the position of a point within the segment.

# Examples

Estimate arc length of segment ``[i, i + 1]``, given by ``ℓ = ∫_{t_{i}}^{t_{i + 1}} |∂_t \bm{X}(t)| \, \mathrm{d}t``:

```julia
quad = GaussLegendre(4)  # quadrature rule
ℓ = integrate(f, i, quad) do f, i, ζ
    norm(f(i, ζ, Derivative(1)))  # = |∂ₜ𝐗|
end
```

Alternatively:
```julia
quad = GaussLegendre(4)  # quadrature rule
s = Segment(f, i)
ℓ = integrate(s, quad) do s, ζ
    norm(s(ζ, Derivative(1)))  # = |∂ₜ𝐗|
end
```
"""
function integrate(
        integrand::F, f::AbstractFilament, i::Int, quad::AbstractQuadrature;
        _args = (f, i),  # integrand arguments; used internally
    ) where {F}
    ζs, ws = quadrature(quad)
    ts = knots(f)
    Δt = ts[i + 1] - ts[i]
    Δt * sum(eachindex(ws)) do j
        @inline
        fx = @inline integrand(_args..., ζs[j])
        ws[j] * fx
    end
end

function integrate(integrand::F, s::Segment, quad::AbstractQuadrature) where {F}
    integrate(integrand, s.f, s.i, quad; _args = (s,))
end

"""
    integrate(integrand::Function, f::AbstractFilament, quad::AbstractQuadrature)

Estimate integral over a whole filament.

The integral is computed as the sum of the integrals over each filament segment.
The given quadrature rule is applied over each segment, meaning that the total number of
evaluations is given by the number of segments multiplied by the length of the quadrature
rule.

The signature of the integrated function must be
`integrand(f::AbstractFilament, i::Int, ζ::Real)`, where `i` is the index of the segment of
interest. See below for some examples.

## Examples

Estimate the total length of a closed filament, ``L = ∮ |∂_t \\bm{X}(t)| \\, \\mathrm{d}t``:

```julia
quad = GaussLegendre(4)  # quadrature rule

integrand(f, i, ζ) = norm(f(i, ζ, Derivative(1)))  # = |∂ₜ𝐗|

# Here `f` is an existent filament.
L = integrate(integrand, f, quad)

# Or, more conveniently, using a `do` block to define an anonymous function.
L = integrate(f, quad) do ff, i, ζ
    norm(ff(i, ζ, Derivative(1)))  # = |∂ₜ𝐗|
end
```
"""
function integrate(integrand::F, f::AbstractFilament, quad::AbstractQuadrature) where {F}
    sum(segments(f)) do seg
        integrate(integrand, seg.f, seg.i, quad)
    end
end

## ====================================================================================== ##
#  No quadrature case: the idea is to only evaluate quantities on filament nodes,
#  avoiding any interpolation.
## ====================================================================================== ##

# This is an internal type which is just used to replace the interpolation syntax
# f(i, ζ, args...) by the on-node evaluation syntax f[i, args...].
# When integrating, we consider that we evaluate in the middle of each segment (at ζ = 1/2).
# Moreover, the location and the derivative at that point are approximated from the
# locations at the two extremities of the segment, according to a straight-segment
# approximation.
struct NoQuadFilament{T, F <: AbstractFilament{T}} <: AbstractFilament{T}
    f :: F
end

@propagate_inbounds (u::NoQuadFilament)(i::Int, ζ::Number) =
    u(i, ζ, Derivative(0))

@propagate_inbounds function (u::NoQuadFilament)(i::Int, ζ::Number, ::Derivative{0})
    (; f,) = u
    # Note: precision could be improved by shifting the "straight" segment in the direction
    # opposite to the curvature vector.
    # See McGreivy et al, Phys. Plasmas 28, 082111 (2021) [https://doi.org/10.1063/5.0058014].
    # However estimating curvatures increases the computation cost, and we want to keep
    # NoQuadrature as cheap as possible.
    (f[i] + f[i + 1]) ./ 2
end

@propagate_inbounds function (u::NoQuadFilament)(i::Int, ζ::Number, ::Derivative{1})
    (; f,) = u
    ts = knots(f)
    (f[i + 1] - f[i]) ./ (ts[i + 1] - ts[i])
end

# Not sure if this is used...
@propagate_inbounds (u::NoQuadFilament)(i, ζ::Number, d) = (u.f[i, d] + u.f[i + 1, d]) ./ 2

function integrate(
        integrand::F, f::AbstractFilament, i::Int, quad::NoQuadrature;
        _args = (NoQuadFilament(f), i),
    ) where {F}
    ts = knots(f)
    args = _noquad_replace_args(_args...)  # replaces filaments by NoQuadFilament
    Δt = ts[i + 1] - ts[i]
    ζ = 0.5  # generally unused
    Δt * integrand(args..., ζ)
end

# These are the only 2 possible options for _args:
_noquad_replace_args(f::NoQuadFilament, i::Int) = (f, i)
_noquad_replace_args(s::Segment) = (Segment(NoQuadFilament(s.f), s.i),)  # replace AbstractFilament by NoQuadFilament
