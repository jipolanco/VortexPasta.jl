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

integrand(f, i, ζ) = f(i, ζ, Derivative(1))  # = |∂ₜ𝐗|

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
# Technically this is only correct when ζ = 0, but this is not checked.
struct NoQuadFilament{F <: AbstractFilament}
    f :: F
end

(u::NoQuadFilament)(i::Int, ζ::Number, args...) = u.f[i, args...]

function integrate(
        integrand::F, f::AbstractFilament, i::Int, quad::NoQuadrature;
        _args = (),
    ) where {F}
    ts = knots(f)
    args = _noquad_replace_args(_args...)  # replaces any AbstractFilament by a NoQuadFilament
    Δt = Quadratures.increment(quad, ts, i)
    ζ = 0  # unused but needed
    Δt * integrand(args..., ζ)
end

_noquad_replace_args(f::AbstractFilament, args...) = (NoQuadFilament(f), _noquad_replace_args(args...)...)
_noquad_replace_args(x::Any, args...) = (x, _noquad_replace_args(args...)...)
_noquad_replace_args() = ()
