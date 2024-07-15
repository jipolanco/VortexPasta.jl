export stretching_rate

@doc raw"""
    stretching_rate(iter::VortexFilamentSolver; quad = nothing) -> Real
    stretching_rate(fs, vs; quad = nothing) -> Real

Compute stretching rate of one or more vortices.

The stretching rate has units ``L T^{-1}`` and is given by:

```math
\frac{\mathrm{d} \mathcal{L}}{\mathrm{d} t}
= ∮ \frac{∂\bm{v}}{∂ξ} ⋅ \mathrm{d}\bm{s}
= - ∮ \bm{v} ⋅ \bm{s}'' \, \mathrm{d}ξ
```

where ``ξ`` is the arc length, ``\bm{v}(ξ)`` the local filament velocity, ``\bm{s}''(ξ)``
the local curvature vector, and $\mathcal{L}$ the instantaneous vortex length.
The last equality is obtained using integration by parts.

In the implementation, the last expression is the one used to compute the stretching rate.

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
    _stretching_rate(quad, f, vs)
end

function _stretching_rate(quad::Nothing, f::ClosedFilament, vs)
    T = number_type(vs)
    x = zero(T)
    ts = Filaments.knots(f)
    for i ∈ eachindex(f, vs)
        v⃗ = vs[i]
        s⃗′ = f[i, Derivative(1)]
        s⃗″ = f[i, Derivative(2)]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        s′³ = sqrt(sum(abs2, s⃗′))^3  # = |s⃗′|³
        # Note: the integrand is (v⃗ ⋅ ρ⃗) * |s⃗′|, where ρ⃗ is the curvature vector and
        # |s⃗′| = dζ/dt is the conversion from arc length ζ to arbitrary
        # parametrisation t.
        x += (v⃗ ⋅ (s⃗′ × (s⃗′ × s⃗″))) * (δt / s′³)
    end
    x
end

# Case with quadratures (requires interpolating the velocity along filaments)
function _stretching_rate(quad, f::ClosedFilament, vs)
    _stretching_rate(isinterpolable(vs), quad, f, vs)
end

function _stretching_rate(::IsInterpolable{true}, quad, f::ClosedFilament, vs)
    T = number_type(vs)
    x = zero(T)
    for i ∈ eachindex(segments(f))
        x += integrate(f, i, quad) do f, i, ζ
            v⃗ = vs(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))
            s⃗″ = f(i, ζ, Derivative(2))
            s′³ = sqrt(sum(abs2, s⃗′))^3  # = |s⃗′|³
            # Note: the integrand is (v⃗ ⋅ ρ⃗) * |s⃗′|, where ρ⃗ is the curvature vector and
            # |s⃗′| = dζ/dt is the conversion from arc length ζ to arbitrary
            # parametrisation t.
            (v⃗ ⋅ (s⃗′ × (s⃗′ × s⃗″))) / s′³
        end
    end
    x
end

function _stretching_rate(::IsInterpolable{false}, quad, f::ClosedFilament, vs)
    T = number_type(vs)
    x = zero(T)

    method = Filaments.discretisation_method(f)
    ts = Filaments.knots(f)
    M = Filaments.npad(method)
    Np = length(f)

    # We use Bumper to avoid allocations managed by Julia's garbage collector.
    buf = Bumper.default_buffer()

    @no_escape buf begin
        V = eltype(vs)
        data = @alloc(V, Np + 2M)
        cs = PaddedVector{M}(data)
        nderiv = Filaments.required_derivatives(method)
        cderiv = ntuple(Val(nderiv)) do _
            local data = @alloc(V, Np + 2M)
            PaddedVector{M}(data)
        end
        # We interpolate the velocity v⃗, which we then integrate along filaments.
        coefs = Filaments.init_coefficients(method, cs, cderiv)
        copyto!(cs, vs)
        Filaments.compute_coefficients!(coefs, ts)
        for i ∈ eachindex(segments(f))
            x += integrate(f, i, quad) do f, i, ζ
                s⃗′ = f(i, ζ, Derivative(1))
                s⃗″ = f(i, ζ, Derivative(2))
                s′³ = sqrt(sum(abs2, s⃗′))^3  # = |s⃗′|³
                v⃗ = Filaments.evaluate(coefs, ts, i, ζ, Derivative(0))
                # Note: the integrand is (v⃗ ⋅ ρ⃗) * |s⃗′|, where ρ⃗ is the curvature vector and
                # |s⃗′| = dζ/dt is the conversion from arc length ζ to arbitrary
                # parametrisation t.
                (v⃗ ⋅ (s⃗′ × (s⃗′ × s⃗″))) / s′³
            end :: T
        end
    end

    x
end
