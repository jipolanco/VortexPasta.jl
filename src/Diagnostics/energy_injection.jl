using LinearAlgebra: ⋅, ×

export energy_injection_rate

@doc raw"""
    energy_injection_rate(iter::VortexFilamentSolver, [vL]; quad = nothing) -> Real
    energy_injection_rate(fs, vL, vs, p::ParamsBiotSavart; quad = nothing) -> Real

Compute energy injection rate from current filament velocities.

This is typically the energy injected by an advecting vortex velocity `vL`, given the
self-induced (Biot–Savart) velocity `vs`.
It can be negative if `vL` actually dissipates energy.
It does not include the contributions to the dissipation by vortex reconnections or
numerical resolution effects.

# Arguments

- `fs`: vortex filament locations;

- `vL`: velocity of filament nodes (including forcing or external velocities);

- `vs`: self-induced filament velocity (from Biot–Savart law);

- `p::ParamsBiotSavart`: Biot–Savart parameters (see [`ParamsBiotSavart`](@ref));

- `quad`: an optional quadrature rule (such as [`GaussLegendre`](@ref VortexPasta.Quadratures.GaussLegendre)) for accurate line integrals.

# Extended help

## Definition

The energy injection rate (in ``L^2 T^{-3}`` units) can be expressed as

```math
ε_{\text{inj}} = \frac{Γ}{V} ∮ \left[ \bm{s}' × \bm{v}_{\text{s}} \right] ⋅ \bm{v}_{\text{L}} \, \mathrm{d}ξ,
```

where ``\bm{v}_{\text{s}}`` is the self-induced vortex velocity (according to the Biot–Savart law)
and ``\bm{v}_{\text{L}}`` is the actual velocity used to advect the vortices. Typically, the
latter can be written as

```math
\bm{v}_{\text{L}} = \bm{v}_{\text{s}} + \bm{v}_{\text{ext}}
```

where ``\bm{v}_{\text{ext}}`` is an externally applied velocity (for example, representing
the interaction with a normal fluid). Note that ``\bm{v}_{\text{s}}`` doesn't contribute
to energy injection, so that energy is conserved (excluding other dissipative effects) in
the absence of an external velocity.
"""
function energy_injection_rate(
        fs::VectorOfFilaments,
        vL::SetOfFilamentsData,
        vs::SetOfFilamentsData,
        p::ParamsBiotSavart{T};
        kws...
    ) where {T <: AbstractFloat}
    @assert T === number_type(fs) === number_type(vL) === number_type(vs)
    ε = zero(T)
    for i in eachindex(fs, vL, vs)
        ε += energy_injection_rate(fs[i], vL[i], vs[i], p; kws...)
    end
    ε
end

function energy_injection_rate(
        f::AbstractFilament, vL::SingleFilamentData, vs::SingleFilamentData, p::ParamsBiotSavart{T};
        quad = nothing,
    ) where {T <: AbstractFloat}
    _energy_injection_rate(quad, f, vL, vs, p)
end

# 1. No quadratures (cheaper)
function _energy_injection_rate(quad::Nothing, f, vL, vs, p)
    prefactor = p.Γ / _domain_volume(p.Ls)
    T = eltype(p)
    ε = zero(T)
    ts = knots(f)
    for i in eachindex(f, vL, vs)
        v⃗_s = vs[i]
        v⃗_L = vL[i]
        s⃗′ = f[i, Derivative(1)]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        ε += ((s⃗′ × v⃗_s) ⋅ v⃗_L) * δt
    end
    prefactor * ε
end

# 2. With quadratures (requires interpolating along filaments)
function _energy_injection_rate(quad, f, vL, vs, args...)
    with_interpolable_fields = isinterpolable(vL) & isinterpolable(vs)
    _energy_injection_rate(with_interpolable_fields, quad, f, vL, vs, args...)
end

function _energy_injection_rate(::IsInterpolable{true}, quad, f, vL, vs, p)
    prefactor = p.Γ / _domain_volume(p.Ls)
    T = eltype(p)
    ε = zero(T)
    for i in eachindex(f, vL, vs)
        ε += integrate(f, i, quad) do f, i, ζ
            v⃗_s = vs(i, ζ)
            v⃗_L = vL(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))
            (s⃗′ × v⃗_s) ⋅ v⃗_L
        end :: T
    end
    prefactor * ε
end

# To be implemented (not sure we need it...)
function _energy_injection_rate(::IsInterpolable{false}, quad, f, vL, vs, p)
    error(
        """
        computation of energy injection rate without interpolable fields not yet implemented.
        Try disabling quadratures (passing `quad = nothing`) or using interpolable filament velocities (e.g. described by Filament objects).
        """
    )
end
