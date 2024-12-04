"""
    SyntheticFields

Provides implementations of synthetic vector fields.

These can be used to represent a "normal" fluid velocity field which influences the motion
of vortex lines.
"""
module SyntheticFields

export
    SyntheticVectorField,
    UniformVectorField,
    FourierSyntheticVectorField,
    FourierBandVectorField

using StaticArrays: SVector
using Random: Random, AbstractRNG
using LinearAlgebra: ⋅, ×

"""
    SyntheticVectorField{T, N}

Abstract type representing a synthetic vector field in ``N`` dimensions.

Here `T <: AbstractFloat` is the type of the returned values when evaluating the field at a
position.

A field can be evaluated using the `f(x⃗)` syntax, where `f` is a `SyntheticVectorField` and
`x⃗` is a physical location, returning an `SVector{N, T}`. Here `x⃗` can be an `N`-element
tuple or `SVector`.
"""
abstract type SyntheticVectorField{T <: AbstractFloat, N} end

"""
    UniformVectorField{T, N} <: SyntheticVectorField{T, N}
    UniformVectorField(u⃗)

Represents a uniform (constant) vector field with no spatial fluctuations.

This can be used for instance to represent a uniform counterflow.

The input `u⃗` can be an `SVector{N, T}` or a tuple.
"""
struct UniformVectorField{T, N} <: SyntheticVectorField{T, N}
    u⃗::SVector{N, T}
end

UniformVectorField(u⃗::Tuple) = UniformVectorField(SVector(u⃗))

(field::UniformVectorField)(x⃗) = field.u⃗

"""
    FourierSyntheticVectorField{T, N} <: SyntheticVectorField{T, N}

Abstract type representing a synthetic vector field implemented in Fourier space.

See also [`SyntheticVectorField`](@ref).
"""
abstract type FourierSyntheticVectorField{T, N} <: SyntheticVectorField{T, N} end

@doc raw"""
    FourierBandVectorField{T, N} <: FourierSyntheticVectorField{T, N}
    FourierBandVectorField(undef, Ls::NTuple; kmin, kmax)

Implements a synthetic vector field in Fourier space.

This type is adapted for vector fields described by a relatively small number of non-zero Fourier
modes. The non-zero modes are within a "band" given by ``k_{\min} ≤ |\bm{k}| ≤ k_{\max}``.

One should initialise the Fourier coefficients of the vector field before performing any evaluations.
For this one can call [`SyntheticFields.init_coefficients!`](@ref) after creating the vector field.
After that, one can evaluate the field as described in [`SyntheticVectorField`](@ref).

# Positional arguments

- `undef`: this is Julia's `undef` variable, used here to explicitly indicate that Fourier
  coefficients are not initialised by this function (but only allocated).

- `Ls`: domain period in each direction. For example, `Ls = (2π, 2π, 2π)` for a ``2π``-periodic cubic domain.

# Keyword arguments

- `kmin`: minimum forcing wavenumber (magnitude);

- `kmax`: maximum forcing wavenumber (magnitude).

These should satisfy `0 ≤ kmin ≤ kmax`.
Moreover, one usually wants `0 < kmin` to ensure that the generated field has zero mean value.
"""
struct FourierBandVectorField{T <: AbstractFloat, N}
    qs :: Vector{NTuple{N, Int}}          # forced wave vectors (length Nf) -- normalised to integer values
    cs :: Vector{SVector{N, Complex{T}}}  # Fourier coefficients of velocity field (length Nf)
    Δks :: NTuple{N, T}                   # wavenumber increment in each direction (Δk = 2π/L)
end

# Determine whether a wavevector k⃗ should be discarded to preserve Hermitian symmetry
# (ensuring real-valued fields in physical space). The coefficient associated to a discarded
# wavevector is redundant as it can be directly obtained from that associated to -k⃗ (its
# conjugate).
#
# In 3D, this will discard wavevectors (kx, ky, kz) satisfying either:
#
#  * kx < 0                       -- can be replaced with (-kx, -ky, -kz) with -kx ≥ 0;
#  * kx == 0 && ky < 0            -- can be replaced with (  0, -ky, -kz) with -ky ≥ 0;
#  * kx == 0 && ky == 0 && kz < 0 -- can be replaced with (  0,   0, -kz) with -kz ≥ 0.
#
@inline discard_wavevector(k⃗) = _discard_wavevector(Tuple(k⃗), Val(+1))

@inline function _discard_wavevector(k⃗::Tuple, _sign::Val{sign}) where {sign}
    kx = k⃗[begin] * sign
    kx < 0 && return true
    kx > 0 && return false
    k⃗_tail = Base.tail(k⃗)
    _discard_wavevector(k⃗_tail, _sign)
end

@inline _discard_wavevector(k⃗::Tuple{}, ::Val) = false

# Returns true if the Fourier coefficient associated to the -k⃗ wavevector is not explicitly
# included and should be obtained from that associated to +k⃗ (as its conjugate).
# We simply check if -k⃗ has been discarded.
@inline has_implicit_conjugate(k⃗) = _discard_wavevector(Tuple(k⃗), Val(-1))

function FourierBandVectorField(::UndefInitializer, Ls::NTuple{N, T}; kmin, kmax) where {T, N}
    kmin = T(kmin)
    kmax = T(kmax)
    0 ≤ kmin ≤ kmax || throw(ArgumentError("expected 0 < kmin ≤ kmax"))
    # Determine forced wave vectors k⃗
    # Notation: "q" represents a normalised wavenumber, q = k / Δk (takes integer values)
    qs_force = NTuple{N, Int}[]  # normalised forced wave vectors (negative indices represent negative k's)
    Δks = T(2π) ./ Ls
    qmax_x = ceil(Int, kmax / Δks[1])
    qs_x = 0:qmax_x  # test wavenumbers along first dimension (Hermitian symmetry: discard modes kx < 0)
    qs_tail = ntuple(Val(N - 1)) do d
        qmax_d = ceil(Int, kmax / Δks[d + 1])
        -qmax_d:qmax_d  # test wavenumbers along dimension d + 1
    end
    qs_iter = Iterators.product(qs_tail...)
    kmin² = kmin * kmin
    kmax² = kmax * kmax
    for q⃗_tail ∈ qs_iter, qx ∈ qs_x
        q⃗ = (qx, q⃗_tail...)
        k⃗ = q⃗ .* Δks
        k² = sum(abs2, k⃗)
        if kmin² ≤ k² ≤ kmax² && !discard_wavevector(q⃗)
            push!(qs_force, q⃗)
        end
    end
    cs = similar(qs_force, SVector{N, Complex{T}})
    FourierBandVectorField(qs_force, cs, Δks)
end

# Generate divergence-free vector (could be generalised to 2D)
function generate_divfree(u⃗::T, k⃗::SVector{3}) where {T <: SVector{3}}
    k² = sum(abs2, k⃗)
    iszero(k²) ? u⃗ : T(u⃗ - ((k⃗ ⋅ u⃗) / k²) * k⃗)
end

"""
    SyntheticFields.init_coefficients!([rng::AbstractRNG,] f::FourierBandVectorField, u_rms::Real)

Initialise Fourier coefficients of vector field.

The "active" Fourier modes are chosen randomly following a normal distribution such that the
resulting velocity components have an rms value given by `u_rms` (on average).
The generated field is divergence-free (in 3D).

# Positional arguments

- `rng` (optional): random number generator for determining Fourier coefficients of the
  forcing velocity field. Default value is `Random.default_rng()`.

- `f`: field to be initialised (of type [`FourierBandVectorField`](@ref)).

- `u_rms`: typical magnitude (rms value) of _each_ vector component.

The total rms value of the vector field is `N * u_rms` where `N` is the number of dimensions.
Note that the rms value of a single component may slightly differ from `u_rms`.
"""
function init_coefficients!(rng::AbstractRNG, f::FourierBandVectorField{T, N}, u_rms::Real) where {T, N}
    (; qs, cs, Δks,) = f
    V = eltype(cs)
    sum2 = zero(T)
    for i ∈ eachindex(qs, cs)
        u⃗ = randn(rng, V)  # all modes have the same std
        q⃗ = qs[i]
        k⃗ = SVector(q⃗ .* Δks)
        u⃗ = generate_divfree(u⃗, k⃗)  # make it divergence-free (orthogonal to k⃗)
        cs[i] = u⃗
        factor = 1 + has_implicit_conjugate(k⃗)  # include -k⃗ in the sum if kx ≠ 0 (Hermitian symmetry)
        sum2 += factor * sum(abs2, u⃗)
    end
    # Ensure that modes (0, ky, kz...) respect Hermitian symmetry (this is not automatic).
    # We want u_rms^2 to be the (average) variance of a *single* velocity component.
    # Therefore, N * u_rms^2 is the variance of the full field.
    normalisation = sqrt(N * T(u_rms)^2 / sum2)
    @. cs = cs * normalisation
    f
end

init_coefficients!(f::FourierBandVectorField, u_rms::Real; kws...) =
    init_coefficients!(Random.default_rng(), f, u_rms; kws...)

# Copy coefficients to regular tuple of arrays (one array = one vector component), which can
# be used for FFTs or other verifications.
# Currently this is only used in tests.
# Note: this may fail or give wrong results if the ûs arrays are not large enough.
function to_fourier_grid!(ûs::NTuple{N, AbstractArray}, field::FourierBandVectorField{T, N}) where {T, N}
    (; qs, cs,) = field
    Ms = size(ûs[1])
    for u ∈ ûs
        fill!(u, 0)
    end
    for i ∈ eachindex(qs, cs)
        q⃗ = qs[i]
        js = map(q⃗, Ms) do q, M
            ifelse(q ≥ 0, q + 1, M + q + 1)  # destination index
        end
        û = cs[i]
        for n ∈ eachindex(û)
            ûs[n][js...] = û[n]
        end
        # Also include implicit conjugate if kx = 0.
        # For example, the output field includes the wavenumber k⃗ = (0, ky, kz) for ky < 0, but
        # this wavenumber is not included in the FourierBandVectorField since it is
        # redundant with -k⃗ (which is included).
        if iszero(q⃗[1]) && has_implicit_conjugate(q⃗)
            js_conj = ntuple(Val(N - 1)) do d
                q = -q⃗[d + 1]
                M = Ms[d + 1]
                ifelse(q ≥ 0, q + 1, M + q + 1)  # destination index
            end
            for n ∈ eachindex(û)
                ûs[n][js[1], js_conj...] = conj(û[n])
            end
        end
    end
    ûs
end

function evaluate_direct(field::FourierBandVectorField{T, N}, x⃗) where {T, N}
    @assert length(x⃗) == N
    (; qs, cs, Δks,) = field
    w⃗ = zero(SVector{N, T})
    for i ∈ eachindex(qs, cs)
        u⃗ = cs[i]
        q⃗ = SVector(qs[i])  # normalised wave vector (integer values)
        k⃗ = q⃗ .* Δks
        θ::T = k⃗ ⋅ x⃗
        s, c = sincos(θ)
        factor = 1 + has_implicit_conjugate(k⃗)  # include -k⃗ in the sum if kx ≠ 0 (Hermitian symmetry)
        w⃗ = muladd(factor, real(u⃗) * c - imag(u⃗) * s, w⃗)
    end
    w⃗
end

# Evaluate field at location
(f::FourierBandVectorField)(x⃗) = evaluate_direct(f, x⃗)

end
