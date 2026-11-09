"""
    kinetic_energy(iter::PointVortexSolver) -> Real
    kinetic_energy(vortices::PointVortexState, ψs::AbstractVector{<:Real}, p::ParamsPointVortex) -> Real

Compute kinetic energy per unit mass of 2D velocity field induced by a set of point vortices.

The kinetic energy is obtained from the streamfunction values `ψs` evaluated at all vortex
positions.
"""
function kinetic_energy(vortices::PointVortexState, ψs::AbstractVector{<:Real}, p::ParamsPointVortex)
    (; Ls, Γ) = p
    (; signs) = vortices
    T = typeof(Γ)
    E = zero(T)
    for i in eachindex(ψs, signs)
        @inbounds E += T(BiotSavart.to_real(signs[i]) * ψs[i])
    end
    E = E * T(Γ / 2 / prod(Ls))
    E
end

@doc raw"""
    vortex_impulse(iter::PointVortexSolver) -> Real
    vortex_impulse(vortices::PointVortexState) -> Real

Estimate normalised impulse of a set of point vortices.

This function returns the 2D normalised impulse defined as

```math
\bm{p} = \sum_j s_j \bm{x}_j × \hat{z}
```

where ``s_j \in \{-1, +1\}`` is the vortex sign.

(see e.g. Weiss & McWilliams, Phys. Fluids **3**, 835, 1991)
"""
function vortex_impulse(vortices::PointVortexState)
    (; points, signs) = vortices
    V = eltype(points)
    p⃗ = zero(V)
    @inbounds for i in eachindex(points, signs)
        x, y = points[i]
        p⃗ += BiotSavart.to_real(signs[i]) * V(y, -x)
    end
    p⃗
end
