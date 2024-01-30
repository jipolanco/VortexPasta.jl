export helicity

using LinearAlgebra: ⋅

@doc raw"""
    helicity(iter::VortexFilamentSolver; quad = nothing)
    helicity(fs, vs, Γ; quad = nothing)

Compute helicity of a vortex configuration.

The [helicity](https://en.wikipedia.org/wiki/Hydrodynamical_helicity) is defined as:

```math
\mathcal{H}
= ∫ \bm{u}(\bm{x}) ⋅ \bm{ω}(\bm{x}) \mathrm{d}^3\bm{x}
= Γ ∮ \bm{v}(\bm{s}) ⋅ \mathrm{d}\bm{s}
```

where ``\bm{u}(\bm{x})`` and ``\bm{ω}(\bm{x}) = \bm{∇} × \bm{u}(\bm{x})`` are the velocity
and vorticity fields, and ``\bm{v}(\bm{s})`` is the velocity of the set of vortex filaments.

# Arguments

- `fs`: list of vortex filaments (or a single `AbstractFilament`);
- `vs`: list of velocities on filament nodes (or a single vector of velocities on a filament);
- `Γ`: vortex circulation.

See e.g. [`kinetic_energy_from_streamfunction`](@ref) for the meaning of the optional `quad`
keyword argument.

"""
function helicity end

function helicity(fs::VectorOfFilaments, vs::SetOfFilamentsData, Γ::Real; kws...)
    T = eltype(eltype(eltype(fs)))
    @assert T <: AbstractFloat
    H = zero(T)
    for (f, v) ∈ zip(fs, vs)
        H += helicity(f, v, Γ; kws...)
    end
    H
end

function helicity(f::AbstractFilament, vs::SingleFilamentData, Γ::Real; quad = nothing)
    _helicity(quad, f, vs, Γ)
end

function _helicity(::Nothing, f::AbstractFilament, vs::SingleFilamentData, Γ::Real)
    T = eltype(eltype(f))
    @assert T <: AbstractFloat
    H = zero(T)
    ts = knots(f)
    for i ∈ eachindex(f, vs)
        v⃗ = vs[i]
        s⃗′ = f[i, Derivative(1)]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        H += (v⃗ ⋅ s⃗′) * δt
    end
    T(Γ) * H
end

# Case with quadratures (requires interpolating the velocity along filaments)
function _helicity(quad, f::AbstractFilament, vs::SingleFilamentData, Γ::Real)
    T = eltype(eltype(f))
    @assert T <: AbstractFloat
    H = zero(T)

    method = Filaments.discretisation_method(f)
    ts = Filaments.knots(f)
    M = Filaments.npad(method)
    Np = length(f)

    # We use Bumper to avoid allocations managed by Julia's garbage collector.
    buf = Bumper.default_buffer()

    @no_escape buf begin
        data = @alloc(T, Np + 2M)
        cs = PaddedVector{M}(data)
        nderiv = Filaments.required_derivatives(method)
        cderiv = ntuple(Val(nderiv)) do _
            local data = @alloc(T, Np + 2M)
            PaddedVector{M}(data)
        end
        # Note: we interpolate the tangent component of the velocity (i.e. v⃗ ⋅ s⃗′),
        # which we then integrate along filaments.
        coefs = Filaments.init_coefficients(method, cs, cderiv)
        for i ∈ eachindex(cs, f, vs)
            @inbounds cs[i] = vs[i] ⋅ f[i, Derivative(1)]
        end
        Filaments.compute_coefficients!(coefs, ts)
        for i ∈ eachindex(segments(f))
            H += integrate(f, i, quad) do f, i, ζ
                Filaments.evaluate(coefs, ts, i, ζ)
            end :: T
        end
    end

    T(Γ) * H
end
