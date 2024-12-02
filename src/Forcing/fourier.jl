export ConstantFourierNormalFluidForcing,
       FourierNormalFluidForcing

"""
    FourierNormalFluidForcing <: NormalFluidForcing

Normal fluid forcing implemented in Fourier space.

See also [`NormalFluidForcing`](@ref).
"""
abstract type FourierNormalFluidForcing <: AbstractForcing end

struct FourierForcingData{T <: AbstractFloat, N}
    ks :: Vector{NTuple{N, Int}}          # forced wave vectors (length Nf)
    cs :: Vector{SVector{N, Complex{T}}}  # Fourier coefficients of velocity field (length Nf)
end

function FourierForcingData{T, N}(; kmin::T, kmax::T) where {T, N}
    0 ≤ kmin ≤ kmax || throw(ArgumentError("expected 0 < kmin ≤ kmax"))
    # Determine forced wave vectors k⃗ (taking integer values)
    ks_force = NTuple{N, Int}[]
    kmax_int = ceil(Int, kmax)
    ks_i = -kmax_int:kmax_int  # test values in each direction
    ks_iter = Iterators.product(ntuple(_ -> ks_i, Val(N - 1))...)  # iterator over all test wavenumbers k⃗ (dimensions 2:N)
    ks_x = 0:kmax_int  # Hermitian symmetry: discard modes kx < 0 (dimension 1)
    kmin² = max(kmin * kmin, eps(kmin))  # make sure the zero mode is not included
    kmax² = kmax * kmax
    for k⃗_tail ∈ ks_iter, kx ∈ ks_x
        k⃗ = (kx, k⃗_tail...)
        k² = sum(abs2, k⃗)
        if kmin² ≤ k² ≤ kmax²
            push!(ks_force, k⃗)
        end
    end
    cs = similar(ks_force, SVector{N, Complex{T}})
    FourierForcingData(ks_force, cs)
end

# Set random Fourier coefficients of the forcing.
function init_coefficients!(rng::AbstractRNG, data::FourierForcingData{T, N}, u_rms::T) where {T, N}
    (; ks, cs,) = data
    V = eltype(cs)
    sum2 = zero(u_rms)
    for i ∈ eachindex(ks, cs)
        u⃗ = randn(rng, V)  # all modes have the same std
        k⃗ = SVector(Tuple(ks[i]))
        u⃗ = u⃗ × k⃗  # make it divergence-free (orthogonal to k⃗)
        cs[i] = u⃗
        factor = 2 - iszero(k⃗[1])  # include -k⃗ in the sum if kx ≠ 0 (Hermitian symmetry)
        @show factor, k⃗
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
    ConstantFourierNormalFluidForcing([rng::AbstractRNG]; α, α′ = 0, vn_rms = 1.0, kmin = 0.5, kmax = 1.5)

Forcing via a normal fluid flow represented in Fourier space.

The normal fluid velocity field is taken as constant. Its Fourier modes are chosen randomly
following a normal distribution such that the resulting velocity components have an rms value
given by `vn_rms` (on average). Moreover, the generated field is divergence-free.

See [`NormalFluidForcing`](@ref) for the general form of this type of forcing.

See also [`FourierNormalFluidForcing`](@ref).

# Positional arguments

- `rng` (optional): random number generator for determining Fourier coefficients of the
  forcing velocity field. Default value is `Random.default_rng()`.

# Keyword arguments

- `α`: Magnus force coefficient;

- `α′ = 0`: drag force coefficient;

- `vn_rms = 1.0`: typical magnitude (rms value) of normal fluid velocity fluctuations;

- `kmin = 0.5`: minimum forcing wavenumber (magnitude);

- `kmax = 1.5`: maximum forcing wavenumber (magnitude).

For simplicity, the minimum/maximum wavenumbers assume a ``(2π)^3``-periodic domain, so that
the wavenumber components ``k_i`` take integer values. The values of `kmin` and `kmax`
actually fix the number of vector wavenumbers ``\bm{k}`` which will be forced, independently
of the actual domain period.

Also, the spectrum of the imposed velocity field is roughly flat in the forcing range.
The magnitude of each Fourier mode is only determined by the value of `vn_rms` (and the
number of forced Fourier modes).
"""
struct ConstantFourierNormalFluidForcing{T <: AbstractFloat} <: FourierNormalFluidForcing
    α  :: T  # Magnus force coefficient
    α′ :: T  # drag force coefficient
    vn_rms :: T
    data :: FourierForcingData{T, 3}  # 3 dimensions only
end

function ConstantFourierNormalFluidForcing(
        rng::AbstractRNG = Random.default_rng();
        α, α′ = 0, vn_rms = 1.0, kmin = 0.5, kmax = 1.5,
    )
    T = float(typeof(α))
    data = FourierForcingData{T, 3}(kmin = T(kmin), kmax = T(kmax))
    init_coefficients!(rng, data, T(vn_rms))
    ConstantFourierNormalFluidForcing(α, T(α′), T(vn_rms), data)
end
