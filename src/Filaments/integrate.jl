"""
    integrate(func::Function, f::AbstractFilament, i::Int, quad::AbstractQuadrature)
    integrate(func::Function, s::Segment, quad::AbstractQuadrature)

Estimate integral along a filament segment using the chosen quadrature.

Integration is performed in the segment `(f[i], f[i + 1])`.

The function `func(ζ)` takes a single argument ``ζ ∈ [0, 1]`` which corresponds
to the position of a point within the segment.

# Examples

Estimate arc length of segment ``[i, i + 1]``, given by ``ℓ = ∫_{t_{i}}^{t_{i + 1}} |∂_t \\bm{X}(t)| \\, \\mathrm{d}t``:

```julia
quad = GaussLegendre(4)  # quadrature rule
ℓ = integrate(f, i, quad) do ζ
    norm(f(i, ζ, Derivative(1)))  # = |∂ₜ𝐗|
end
```
"""
function integrate(
        func::F, f::AbstractFilament, i::Int, quad::AbstractQuadrature;
        _args = _prepended_args(func, f, i),  # used internally
    ) where {F}
    ζs, ws = quadrature(quad)
    ts = knots(f)
    Δt = ts[i + 1] - ts[i]
    Δt * sum(eachindex(ws)) do j
        @inline
        fx = @inline func(_args..., ζs[j])
        ws[j] * fx
    end
end

# Depending on the number of arguments required by the function, prepend some arguments.
function _prepended_args(func::F, f::AbstractFilament, i) where {F}
    if hasmethod(func, Tuple{Any,Any,Any})  # func(f, i, ζ)
        (f, i)
    elseif hasmethod(func, Tuple{Any})      # func(ζ)
        ()
    end
end

function _prepended_args(func::F, s::Segment) where {F}
    if hasmethod(func, Tuple{Any,Any})  # func(seg, ζ)
        (s,)
    else
        _prepended_args(func, s.f, s.i)
    end
end

function integrate(func::F, s::Segment, quad::AbstractQuadrature) where {F}
    (; f, i,) = s
    integrate(func, f, i, quad; _args = _prepended_args(func, s))
end

"""
    integrate(func::Function, f::AbstractFilament, quad::AbstractQuadrature)

Estimate integral over a whole filament.

The integral is computed as the sum of the integrals over each filament segment.
The given quadrature rule is applied over each segment, meaning that the total number of
evaluations is given by the number of segments multiplied by the length of the quadrature
rule.

The integrated function can take either 2 or 3 arguments, in which case they are
respectively interpreted as `func(seg::Segment, ζ::Real)` or
`func(f::AbstractFilament, i::Int, ζ::Real)`.

## Examples

Estimate the total length of a closed filament, ``L = ∮ |∂_t \\bm{X}(t)| \\, \\mathrm{d}t``:

```julia
quad = GaussLegendre(4)  # quadrature rule

# Alternative 1
L = integrate(f, quad) do seg, ζ
    norm(seg(ζ, Derivative(1)))  # = |∂ₜ𝐗|
end

# Alternative 2
L = integrate(f, quad) do ff, ii, ζ
    norm(ff(ii, ζ, Derivative(1)))  # = |∂ₜ𝐗|
end
```
"""
function integrate(func::F, f::AbstractFilament, quad::AbstractQuadrature) where {F}
    sum(segments(f)) do seg
        integrate(func, seg, quad)
    end
end
