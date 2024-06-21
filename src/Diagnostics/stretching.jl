export stretching_rate

@doc raw"""
    stretching_rate(fs, vs; quad = nothing) -> Real

Compute stretching rate of one or more vortices.

The stretching rate has units ``L T^{-1}`` and is given by

```math
\frac{\mathrm{d} \mathcal{L}}{\mathrm{d} t}
= ∮ \frac{∂\bm{v}}{∂ξ} ⋅ \mathrm{d}\bm{s}
```

where ``ξ`` is the arc length, ``\bm{v}(ξ)`` the local filament velocity, and $\mathcal{L}$
the instantaneous vortex length.

The velocity derivatives are obtained using the same interpolation method used to represent
the filaments `fs`.

## Mandatory arguments

- `vs`: velocity values at filament nodes;

- `fs`: vortex filament locations.

## Optional keyword arguments

- `quad = nothing`: optional quadrature rule (e.g. `quad = GaussLegendre(4)`) used to
  evaluate line integrals. If `nothing`, only values at nodes are used (cheaper).
"""
function stretching_rate end

function stretching_rate(fs::VectorOfFilaments, vs::SetOfFilamentsData; kws...)
    T = number_type(vs)
    x = zero(T)
    for i ∈ eachindex(fs, vs)
        x += stretching_rate(fs[i], vs[i]; kws...)
    end
    x
end

function stretching_rate(f::AbstractFilament, vs::SingleFilamentData; quad = nothing)
    # We interpolate the filament velocities using the same knots `ts` of the filaments.
    # This allows us to obtain the spatial derivative of the velocity `v⃗'(τ)` along the filament.
    method = Filaments.discretisation_method(f)
    ts = Filaments.knots(f)
    M = Filaments.npad(method)
    Np = length(f)

    # We use Bumper to avoid allocations managed by Julia's garbage collector.
    buf = Bumper.default_buffer()

    local x

    @no_escape buf begin
        V = eltype(vs)
        @assert V <: Vec3
        data = @alloc(V, Np + 2M)
        cs = PaddedVector{M}(data)
        nderiv = max(1, Filaments.required_derivatives(method))  # we need at least 1 derivative
        cderiv = ntuple(Val(nderiv)) do _
            local data = @alloc(V, Np + 2M)
            PaddedVector{M}(data)
        end
        # We interpolate the velocity v⃗, which we then derive and integrate along filaments.
        coefs = Filaments.init_coefficients(method, cs, cderiv)
        copyto!(cs, vs)
        Filaments.compute_coefficients!(coefs, ts)
        x = _stretching_rate(quad, f, coefs, ts)
    end

    x
end

function _stretching_rate(quad::Nothing, f::ClosedFilament, coefs, ts)
    T = eltype(eltype(coefs))
    @assert T <: AbstractFloat
    x = zero(T)
    for i ∈ eachindex(f)
        # TODO: avoid interpolations if evaluating right on a node? (makes more sense for Hermite interpolations)
        ζ = zero(T)
        v⃗′ = Filaments.evaluate(coefs, ts, i, ζ, Derivative(1))  # spatial derivative (along the filament) of the filament velocity
        t̂ = f[i, UnitTangent()]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        x += (v⃗′ ⋅ t̂) * δt
    end
    x
end

function _stretching_rate(quad, f, coefs, ts)
    T = eltype(eltype(coefs))
    @assert T <: AbstractFloat
    x = zero(T)
    for i ∈ eachindex(segments(f))
        x += integrate(f, i, quad) do f, i, ζ
            v⃗′ = Filaments.evaluate(coefs, ts, i, ζ, Derivative(1))
            t̂ = f(i, ζ, UnitTangent())
            v⃗′ ⋅ t̂
        end :: T
    end
    x
end
