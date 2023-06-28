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
function integrate(func::F, f::AbstractFilament, i::Int, quad::AbstractQuadrature) where {F}
    ζs, ws = quadrature(quad)
    ts = knots(f)
    Δt = ts[i + 1] - ts[i]
    Δt * sum(eachindex(ws)) do j
        @inline
        fx = @inline func(ζs[j])
        ws[j] * fx
    end
end

function integrate(func::F, s::Segment, quad::AbstractQuadrature) where {F}
    (; f, i,) = s
    integrate(func, f, i, quad)
end
