using HCubature: HCubature

export kinetic_energy_from_streamfunction, kinetic_energy_nonperiodic, kinetic_energy

"""
    kinetic_energy(iter::VortexFilamentSolver; quad = nothing)

Compute kinetic energy of velocity field induced by a set of vortex filaments.

This function calls either [`kinetic_energy_from_streamfunction`](@ref) or
[`kinetic_energy_nonperiodic`](@ref) depending on whether the simulation domain is periodic
or not.
"""
function kinetic_energy end

# ======================================================================================== #

# Periodic case

@doc raw"""
    kinetic_energy_from_streamfunction(ψs, fs, Γ, Ls; quad = nothing)
    kinetic_energy_from_streamfunction(iter::VortexFilamentSolver; quad = nothing)

Compute kinetic energy per unit mass (units ``L^2 T^{-2}``) from streamfunction values at
filament nodes in a periodic domain.

In a periodic domain, the kinetic energy per unit mass of the velocity field induced by a
set of vortex filaments can be obtained as:

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

## Non-periodic domains

If the domain is not periodic, that is, if one or more values in `Ls` is [`Infinity`](@ref), then
the domain volume ``V`` in the expression above will be set to 1.
This corresponds to computing an energy *per unit density* (units ``L^5 T^{-2}``) instead of
per unit mass, meaning that one should multiply by the fluid density ``ρ`` (``M L^{-3}``) to
get an actual kinetic energy (``M L^2 T^{-2}``).

However, note that the above definition of the kinetic energy is actually not valid in
non-periodic domains (since it comes from an integration by parts, and the boundary term
can no longer be neglected in this case).
Therefore, in this case one should instead use [`kinetic_energy_nonperiodic`](@ref), which
uses a different definition commonly used for open domains.

"""
function kinetic_energy_from_streamfunction(ψs::AbstractVector, args...; quad = nothing)
    _kinetic_energy_from_streamfunction(quad, ψs, args...)
end

# Case of a set of filaments
function _kinetic_energy_from_streamfunction(
        quad, ψs::SetOfFilamentsData, fs::VectorOfFilaments, Γ, args...,
    )
    T = float(typeof(Γ))
    E = zero(T)
    for i ∈ eachindex(fs, ψs)
        E += _kinetic_energy_from_streamfunction(quad, ψs[i], fs[i], Γ, args...)
    end
    E
end

function _domain_volume(Ls)
    V = prod(Ls)
    if V === Infinity()
        true  # set volume to 1 for infinite domain
    else
        V
    end
end

# Case of a single filament
# 1. No quadratures (cheaper)
function _kinetic_energy_from_streamfunction(
        ::Nothing, ψf::SingleFilamentData, f::AbstractFilament, Γ, Ls,
    )
    prefactor = Γ / (2 * _domain_volume(Ls))
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
    prefactor = Γ / (2 * _domain_volume(Ls))
    E = zero(prefactor)
    # Create Filament object where each node is a streamfunction value instead of a position.
    # This allows to interpolate the streamfunction in-between nodes.
    Xoff = Filaments.end_to_end_offset(f)
    ψ_int = Filaments.change_offset(similar(f), zero(Xoff))
    copy!(Filaments.nodes(ψ_int), ψf)
    Filaments.update_coefficients!(ψ_int; knots = knots(f))
    for i ∈ eachindex(segments(f))
        E += integrate(f, i, quad) do f, i, ζ
            ψ⃗ = ψ_int(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))
            ψ⃗ ⋅ s⃗′
        end
    end
    prefactor * E
end

# Estimate kinetic energy of periodic velocity field in physical space.
# The velocity field is expected to be described by a smooth spatial function.
function kinetic_energy_of_periodic_velocity_field(
        vext::F, Ls;
        rtol = sqrt(eps(eltype(Ls))),
    ) where {F <: Function}
    @assert length(Ls) == 3 && !isinf(prod(Ls))
    bs = Vec3(Ls)
    as = zero(bs)
    f(x⃗) = sum(abs2, vext(x⃗))
    E, err = HCubature.hcubature(f, as, bs; rtol)
    E /= 2 * _domain_volume(Ls)
    E
end

# ======================================================================================== #

# Non-periodic case: use standard energy definition which assumes that the velocity tends to
# 0 at infinity.

@doc raw"""
    kinetic_energy_nonperiodic(vs, fs, Γ; quad = nothing)
    kinetic_energy_nonperiodic(iter::VortexFilamentSolver; quad = nothing)

Compute kinetic energy per unit density (units ``L^5 T^{-2}``) fro velocity values at
filament nodes in an open (non-periodic) domain.

In an open domain, the kinetic energy of the velocity field induced by a set of vortex
filaments can be obtained as:

```math
E = ρ Γ ∮ \bm{v} ⋅ (\bm{s} × \mathrm{d}\bm{s})
```

where ``Γ`` is the vortex circulation and ``ρ`` is the fluid density.
This definition assumes that the velocity field tends to zero far from the vortices (which
is true in open domains, but not in periodic ones).

This function returns the energy *per unit density* ``E / ρ``.

## Mandatory arguments

- `vs`: velocity values at filament nodes;

- `fs`: vortex filament locations;

- `Γ::Real`: quantum of circulation.

Alternatively, one may simply pass a `VortexFilamentSolver`, which contains the current
state of a vortex filament simulation.

## Optional keyword arguments

See [`kinetic_energy_from_streamfunction`](@ref) for details.

## Periodic domains

In fully periodic domains, the [`kinetic_energy_from_streamfunction`](@ref) function should
be used instead to compute the kinetic energy.

"""
function kinetic_energy_nonperiodic(vs::AbstractVector, args...; quad = nothing)
    _kinetic_energy_nonperiodic(quad, vs, args...)
end

function _kinetic_energy_nonperiodic(
        quad, vs::SetOfFilamentsData, fs::VectorOfFilaments, Γ::AbstractFloat,
    )
    E = zero(typeof(Γ))
    for i ∈ eachindex(vs, fs)
        E += _kinetic_energy_nonperiodic(quad, vs[i], fs[i], Γ)
    end
    E
end

# Case without quadratures
function _kinetic_energy_nonperiodic(
        ::Nothing, vs::SingleFilamentData, f::AbstractFilament, Γ,
    )
    E = zero(typeof(Γ))
    ts = knots(f)
    for i ∈ eachindex(f, vs)
        s⃗ = f[i]
        s⃗′ = f[i, Derivative(1)]
        v⃗ = vs[i]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        E += (v⃗ ⋅ (s⃗ × s⃗′)) * δt
    end
    Γ * E
end

# Case with quadratures
function _kinetic_energy_nonperiodic(
        quad, vs::SingleFilamentData, f::AbstractFilament, Γ,
    )
    E = zero(typeof(Γ))
    # Create Filament object where each node is a velocity instead of a position.
    # This allows to interpolate velocities in-between nodes.
    Xoff = Filaments.end_to_end_offset(f)
    v_int = Filaments.change_offset(similar(f), zero(Xoff))
    copy!(Filaments.nodes(v_int), vs)
    Filaments.update_coefficients!(v_int; knots = knots(f))
    for i ∈ eachindex(segments(f))
        E += integrate(f, i, quad) do f, i, ζ
            s⃗ = f(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))
            v⃗ = v_int(i, ζ)
            v⃗ ⋅ (s⃗ × s⃗′)
        end
    end
    Γ * E
end
