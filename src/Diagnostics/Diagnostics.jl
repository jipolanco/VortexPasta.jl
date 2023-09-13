"""
    Diagnostics

Contains tools for computing different diagnostics (total energy, energy spectra, ...) from
simulation data.
"""
module Diagnostics

export kinetic_energy_from_streamfunction

using ..Filaments: Filaments, Derivative, Vec3, knots, segments, integrate
using LinearAlgebra: ⋅

const SingleFilamentData = AbstractVector{<:Vec3}
const SetOfFilamentsData = AbstractVector{<:SingleFilamentData}

@doc raw"""
    kinetic_energy_from_streamfunction(ψs, fs, Γ, Ls; quad = nothing)
    kinetic_energy_from_streamfunction(iter::VortexFilamentSolver; quad = nothing)

Compute kinetic energy per unit mass (units ``L²T⁻²``) from streamfunction values at
filament nodes in a periodic domain.

The kinetic energy per unit mass associated to a set of vortex filaments is defined as:

```math
E = \frac{Γ}{2V} ∮ \bm{ψ}(\bm{s}) ⋅ \mathrm{d}\bm{s}
```

where ``Γ`` is the vortex circulation and ``V`` is the volume of a periodic cell.

## Mandatory arguments

- `ψs`: streamfunction values at filament nodes;

- `fs`: vortex filament locations;

- `Γ::Real`: quantum of circulation;

- `Ls::Tuple = (Lx, Ly, Lz)`: domain period in each direction.

Alternatively, one may simply pass a `VortexFilamentSolver`, which contains the current
state of a vortex filament simulation.

## Optional keyword arguments

- `quad = nothing`: optional quadrature rule (e.g. `quad = GaussLegendre(4)`) used to
  evaluate line integrals. If `nothing`, only values at nodes are used (cheaper). Otherwise,
  if a quadrature rule is passed, interpolations are performed and extra allocations are
  needed.
"""
function kinetic_energy_from_streamfunction(ψs::AbstractVector, args...; quad = nothing)
    _kinetic_energy_from_streamfunction(quad, ψs, args...)
end

function kinetic_energy_from_streamfunction(iter; kws...)
    (; ψs, fs,) = iter
    (; Γ, Ls,) = iter.prob.p.common
    kinetic_energy_from_streamfunction(ψs, fs, Γ, Ls; kws...)
end

# Case of a set of filaments
function _kinetic_energy_from_streamfunction(
        quad, ψs::SetOfFilamentsData, fs::SetOfFilamentsData, Γ, args...,
    )
    sum(eachindex(fs, ψs)) do i
        _kinetic_energy_from_streamfunction(quad, ψs[i], fs[i], Γ, args...)
    end
end

# Case of a single filament
# 1. No quadratures (cheaper)
function _kinetic_energy_from_streamfunction(
        ::Nothing, ψf::SingleFilamentData, f::SingleFilamentData, Γ, Ls,
    )
    prefactor = Γ / (2 * prod(Ls))
    E = zero(prefactor)
    ts = knots(f)
    for i ∈ eachindex(f, ψf)
        ψ⃗ = ψf[i]
        s⃗′ = f[i, Derivative(1)]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        E += (ψ⃗ ⋅ s⃗′) * δt
    end
    prefactor * E
end

# With quadratures (requires interpolation + allocations)
function _kinetic_energy_from_streamfunction(
        quad, ψf::SingleFilamentData, f::SingleFilamentData, Γ, Ls,
    )
    prefactor = Γ / (2 * prod(Ls))
    E = zero(prefactor)
    Xoff = Filaments.end_to_end_offset(f)
    ψ_int = Filaments.change_offset(similar(f), zero(Xoff))
    copy!(Filaments.nodes(ψ_int), ψf)
    Filaments.update_coefficients!(ψ_int; knots = knots(f))
    for i ∈ eachindex(segments(f))
        E += integrate(f, i, quad) do ζ
            ψ⃗ = ψ_int(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))
            ψ⃗ ⋅ s⃗′
        end
    end
    prefactor * E
end

end
