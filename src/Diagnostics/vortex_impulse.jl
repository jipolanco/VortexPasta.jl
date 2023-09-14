export vortex_impulse

@doc raw"""
    vortex_impulse(f; quad = nothing) -> Vec3
    vortex_impulse(iter::VortexFilamentSolver; quad = nothing) -> Vec3

Estimate normalised impulse of one or more vortex filaments.

The vortex impulse is defined as

```math
\bm{p} = \frac{1}{2} ρ Γ ∮ \bm{s} × \mathrm{d}\bm{s}
```

where ``ρ`` is the fluid density and ``Γ`` the circulation about the vortex.
Note that this function returns the *normalised* impulse ``\bm{p} / ρΓ``.
The returned impulse has units of ``L^2`` (an area).

Note that, for a circular vortex ring of radius ``R``, its impulse is ``\bm{p} = ρΓA``
where ``A = π R^2`` is the area enclosed by the ring (and the orientation is equal
to its direction of propagation, i.e. normal to the plane where the ring lives).
"""
function vortex_impulse end

function vortex_impulse(fs::VectorOfFilaments; kws...)
    T = eltype(eltype(fs))
    @assert T <: Vec3
    p⃗ = zero(T)
    for f ∈ fs
        p⃗ = p⃗ + vortex_impulse(f; kws...)
    end
    p⃗
end

vortex_impulse(f::AbstractFilament; quad = nothing) = _vortex_impulse(quad, f)

function _vortex_impulse(quad, f::AbstractFilament)
    T = eltype(f)
    p⃗ = zero(T)
    for seg ∈ segments(f)
        p⃗ = p⃗ + integrate(seg, quad) do ζ
            s⃗ = seg(ζ)
            s⃗′ = seg(ζ, Derivative(1))
            s⃗ × s⃗′
        end
    end
    p⃗ ./ 2
end

function _vortex_impulse(::Nothing, f::AbstractFilament)
    ts = knots(f)
    sum(eachindex(f)) do i
        s⃗ = f[i]
        s⃗′ = f[i, Derivative(1)]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        (s⃗ × s⃗′) * δt
    end ./ 2
end
