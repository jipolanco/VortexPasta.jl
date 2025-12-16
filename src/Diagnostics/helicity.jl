export helicity

using LinearAlgebra: ⋅

@doc raw"""
    helicity(iter::VortexFilamentSolver; quad = nothing) -> Real
    helicity(fs, vs, Γ::Real; quad = nothing) -> Real
    helicity(fs, vs, p::ParamsBiotSavart; quad = nothing) -> Real

Compute helicity of a vortex configuration.

The [helicity](https://en.wikipedia.org/wiki/Hydrodynamical_helicity) is defined as:

```math
\mathcal{H}
= ∫ \bm{u}(\bm{x}) ⋅ \bm{ω}(\bm{x}) \, \mathrm{d}^3\bm{x}
= Γ ∮ \bm{v}(\bm{s}) ⋅ \mathrm{d}\bm{s}
```

where ``\bm{u}(\bm{x})`` and ``\bm{ω}(\bm{x}) = \bm{∇} × \bm{u}(\bm{x})`` are the velocity
and vorticity fields, and ``\bm{v}(\bm{s})`` is the velocity of the set of vortex filaments.

Note that the returned helicity has units of a squared circulation (``Γ^2 ∼ L^4 T^{-2}``).

# Arguments

- `fs`: list of vortex filaments (or a single `AbstractFilament`);

- `vs`: list of velocities on filament nodes (or a single vector of velocities on a filament);

- `p`: Biot–Savart parameters; or

- `Γ`: vortex circulation.

See e.g. [`kinetic_energy_from_streamfunction`](@ref) for the meaning of the optional `quad`
keyword argument.

"""
function helicity end

function helicity(fs::VectorOfFilaments, vs::SetOfFilamentsData, p; quad = nothing, nthreads = Threads.nthreads())  # p could be Γ::Real or a ParamsBiotSavart
    maybe_parallelise_sum(fs, nthreads) do i, inds
        helicity(fs[i], vs[i], p; quad, inds)
    end
end

function helicity(f::AbstractFilament, vs::SingleFilamentData, Γ::Real; quad = nothing, inds = eachindex(f))
    _helicity(quad, f, vs, inds, Γ)
end

function helicity(f::AbstractFilament, vs::SingleFilamentData, p::ParamsBiotSavart; kws...)
    helicity(f, vs, p.Γ; kws...)
end

function _helicity(::Nothing, f::AbstractFilament, vs::SingleFilamentData, inds, Γ::Real)
    @assert eachindex(f) == eachindex(vs)
    checkbounds(f, inds)
    T = number_type(f)
    @assert T <: AbstractFloat
    H = zero(T)
    ts = knots(f)
    @inbounds for i ∈ inds
        v⃗ = vs[i]
        s⃗′ = f[i, Derivative(1)]
        δt = (ts[i + 1] - ts[i - 1]) / 2
        H += (v⃗ ⋅ s⃗′) * δt
    end
    T(Γ) * H
end

# Case with quadratures (requires interpolating the velocity along filaments)
function _helicity(quad, f::AbstractFilament, vs::SingleFilamentData, inds, Γ::Real)
    _helicity(isinterpolable(vs), quad, f, vs, inds, Γ)
end

function _helicity(::IsInterpolable{true}, quad, f::AbstractFilament, vs::SingleFilamentData, inds, Γ::Real)
    @assert eachindex(f) == eachindex(vs)
    checkbounds(f, inds)
    T = number_type(f)
    @assert T <: AbstractFloat
    H = zero(T)
    for i ∈ inds
        H += integrate(f, i, quad) do f, i, ζ
            v⃗ = vs(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))
            v⃗ ⋅ s⃗′
        end :: T
    end
    T(Γ) * H
end

function _helicity(::IsInterpolable{false}, quad, f::AbstractFilament, vs::SingleFilamentData, inds, Γ::Real)
    @assert eachindex(f) == eachindex(vs)
    checkbounds(f, inds)
    T = number_type(f)
    @assert T <: AbstractFloat
    H = zero(T)

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
        # We interpolate the the velocity vector v⃗.
        coefs = Filaments.init_coefficients(method, cs, cderiv)
        copyto!(cs, vs)
        Filaments.compute_coefficients!(coefs, ts)
        for i ∈ inds
            H += integrate(f, i, quad) do f, i, ζ
                v⃗ = Filaments.evaluate(coefs, ts, i, ζ)
                s⃗′ = f(i, ζ, Derivative(1))
                v⃗ ⋅ s⃗′
            end :: T
        end
    end

    T(Γ) * H
end
