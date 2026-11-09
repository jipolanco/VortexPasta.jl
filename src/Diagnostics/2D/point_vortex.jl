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
        @inbounds E += T(signs[i] * ψs[i])
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

where ``s_j ∈ \{-1, +1\}`` is the vortex sign.

(see e.g. Weiss & McWilliams, Phys. Fluids **3**, 835, 1991)
"""
function vortex_impulse(vortices::PointVortexState)
    (; points, signs) = vortices
    V = eltype(points)
    p⃗ = zero(V)
    @inbounds for i in eachindex(points, signs)
        x, y = points[i]
        p⃗ += signs[i] * V(y, -x)
    end
    p⃗
end

function energy_flux(
        cache::LongRangeCache, vortices::PointVortexState, velocities::NamedTuple,
        ks::AbstractVector, params::ParamsPointVortex;
        vs_buf = similar(first(velocities).field),
    )
    (; wavenumbers, state,) = BiotSavart.get_longrange_field_fourier(cache)
    (; Γ, Ls) = params
    (; signs) = vortices
    prefactor = Γ / prod(Ls)
    @assert state.quantity == :velocity
    @assert state.smoothing_scale == 0

    # Initialise flux vectors.
    fluxes = map(_ -> similar(ks), velocities)::NamedTuple

    # Iterate over target wavenumbers.
    vs_filtered = vs_buf
    V = eltype(vs_filtered)  # = Vec2{T}
    for (i, k_flux) in pairs(ks)
        callback = DropLargeWavenumbers(k_flux, wavenumbers)
        output = cache.outputs.velocity
        BiotSavart.interpolate_to_physical!(callback, output, cache)
        BiotSavart.copy_output_values_on_nodes!(vs_filtered, output, cache.pointdata.buf_host)

        # Compute dissipation of large-scale energy by each input velocity
        for n in eachindex(values(velocities), values(fluxes))
            vs = vs_filtered
            vL = velocities[n].field
            s = velocities[n].sign  # -1 if dissipation, +1 if injection
            fluxes[n][i] = s * prefactor * sum(eachindex(vortices)) do j
                @inbounds w⃗ = V(-vs[j][2], +vs[j][1])  # = ẑ × vs
                @inbounds signs[j] * (w⃗ ⋅ vL[j])
            end
        end
    end

    ks, fluxes
end
