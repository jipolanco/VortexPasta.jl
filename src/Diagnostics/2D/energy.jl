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
