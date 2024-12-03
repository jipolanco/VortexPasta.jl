export ConstantFourierNormalFluidForcing,
       FourierNormalFluidForcing

"""
    FourierNormalFluidForcing <: NormalFluidForcing

Normal fluid forcing implemented in Fourier space.

See also [`NormalFluidForcing`](@ref).
"""
abstract type FourierNormalFluidForcing <: AbstractForcing end

struct FourierForcingData{T <: AbstractFloat, N}
    qs :: Vector{NTuple{N, Int}}          # forced wave vectors (length Nf) -- normalised to integer values
    cs :: Vector{SVector{N, Complex{T}}}  # Fourier coefficients of velocity field (length Nf)
    Δks :: NTuple{N, T}                   # wavenumber increment in each direction (Δk = 2π/L)
end

function FourierForcingData(Ls::NTuple{N, T}; kmin::T, kmax::T) where {T, N}
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
    kmin² = max(kmin * kmin, eps(kmin))  # make sure the zero mode is not included
    kmax² = kmax * kmax
    for q⃗_tail ∈ qs_iter, qx ∈ qs_x
        q⃗ = (qx, q⃗_tail...)
        k⃗ = q⃗ .* Δks
        k² = sum(abs2, k⃗)
        if kmin² ≤ k² ≤ kmax²
            push!(qs_force, q⃗)
        end
    end
    cs = similar(qs_force, SVector{N, Complex{T}})
    FourierForcingData(qs_force, cs, Δks)
end

# Set random Fourier coefficients of the forcing.
function init_coefficients!(rng::AbstractRNG, data::FourierForcingData{T, N}, u_rms::T) where {T, N}
    (; qs, cs, Δks,) = data
    V = eltype(cs)
    sum2 = zero(u_rms)
    for i ∈ eachindex(qs, cs)
        u⃗ = randn(rng, V)  # all modes have the same std
        k⃗ = SVector(qs[i] .* Δks)
        u⃗ = u⃗ × k⃗  # make it divergence-free (orthogonal to k⃗)
        cs[i] = u⃗
        factor = 2 - iszero(k⃗[1])  # include -k⃗ in the sum if kx ≠ 0 (Hermitian symmetry)
        sum2 += factor * sum(abs2, u⃗)
    end
    # We want u_rms^2 to be the (average) variance of a *single* velocity component.
    # Therefore, N * u_rms^2 is the variance of the full field.
    normalisation = sqrt(N * u_rms^2 / sum2)
    @. cs = cs * normalisation
    data
end

@doc raw"""
    ConstantFourierNormalFluidForcing <: FourierNormalFluidForcing
    ConstantFourierNormalFluidForcing([rng::AbstractRNG], Ls::NTuple; α, α′ = 0, v_rms, kmin, kmax)

Forcing via a normal fluid flow represented in Fourier space.

The normal fluid velocity field is taken as constant. Its Fourier modes are chosen randomly
following a normal distribution such that the resulting velocity components have an rms value
given by `v_rms` (on average). Moreover, the generated field is divergence-free.

See [`NormalFluidForcing`](@ref) for the general form of this type of forcing.

See also [`FourierNormalFluidForcing`](@ref).

# Positional arguments

- `rng` (optional): random number generator for determining Fourier coefficients of the
  forcing velocity field. Default value is `Random.default_rng()`.

- `Ls`: domain period in each direction. For example, `Ls = (2π, 2π, 2π)` for a ``2π``-periodic cubic domain.

# Keyword arguments

- `α`: Magnus force coefficient;

- `α′ = 0`: drag force coefficient;

- `v_rms`: typical magnitude (rms value) of normal fluid velocity fluctuations;

- `kmin`: minimum forcing wavenumber (magnitude);

- `kmax`: maximum forcing wavenumber (magnitude).

For now, the spectrum of the generated velocity field is roughly flat in the forcing range.
The magnitude of each Fourier mode is only determined by the value of `v_rms` (and the
number of forced Fourier modes).
"""
struct ConstantFourierNormalFluidForcing{T <: AbstractFloat, N} <: FourierNormalFluidForcing
    α  :: T  # Magnus force coefficient
    α′ :: T  # drag force coefficient
    v_rms :: T
    data :: FourierForcingData{T, N}
end

function ConstantFourierNormalFluidForcing(
        rng::AbstractRNG, Ls::NTuple{N, T};
        α, α′ = 0, v_rms, kmin, kmax,
    ) where {N, T <: AbstractFloat}
    data = FourierForcingData(Ls; kmin = T(kmin), kmax = T(kmax))
    init_coefficients!(rng, data, T(v_rms))
    ConstantFourierNormalFluidForcing(T(α), T(α′), T(v_rms), data)
end

ConstantFourierNormalFluidForcing(Ls::NTuple; kws...) =
    ConstantFourierNormalFluidForcing(Random.default_rng(), Ls; kws...)
