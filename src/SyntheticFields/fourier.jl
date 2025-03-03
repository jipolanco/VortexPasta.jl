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

Moreover the state can be saved to an HDF5 file with
[`SyntheticFields.save_coefficients`](@ref) and loaded back with
[`SyntheticFields.load_coefficients!`](@ref).

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
struct FourierBandVectorField{
        T <: AbstractFloat, N,
        NormalisedWaveVectors <: AbstractVector{NTuple{N, Int}},
        FourierCoefficients <: AbstractVector{SVector{N, Complex{T}}},
    } <: FourierSyntheticVectorField{T, N}
    qs::NormalisedWaveVectors  # forced wave vectors (length Nf) -- normalised to integer values
    cs::FourierCoefficients    # Fourier coefficients of velocity field (length Nf)
    Δks::NTuple{N, T}          # wavenumber increment in each direction (Δk = 2π/L)
end

function Base.similar(field::FourierBandVectorField)
    # Note: we copy the wavenumbers; in general we want them to be the same as for the input.
    FourierBandVectorField(copy(field.qs), similar(field.cs), field.Δks)
end

# This allows conversion between CPU and GPU arrays using adapt(...).
@inline function Adapt.adapt_structure(to, field::FourierBandVectorField)
    FourierBandVectorField(
        adapt(to, field.qs),
        adapt(to, field.cs),
        Δks,
    )
end

# Returns the number of elements (vector wavenumbers) describing the field.
# It assumes that qs and cs have the same length.
Base.length(field::FourierBandVectorField) = length(field.cs)

# This is mainly used for CPU <-> GPU transfers.
function Base.copyto!(dst::FourierBandVectorField, src::FourierBandVectorField)
    length(dst) == length(src) || throw(DimensionMismatch("the destination field has the wrong size"))
    copyto!(dst.cs, src.cs)  # copy coefficients only; assume wavenumbers are the same
    nothing
end

function Base.show(io::IO, field::FourierBandVectorField{T, N}) where {T, N}
    (; qs, Δks,) = field
    # Get actual range of wavevector magnitudes
    kext² = extrema(qs) do q⃗
        k⃗ = q⃗ .* Δks
        sum(abs2, k⃗)
    end
    kmin, kmax = round.(sqrt.(kext²); sigdigits = 5)
    Nk = length(qs)
    print(io, "FourierBandVectorField{$T, $N} with $Nk independent Fourier coefficients in |k⃗| ∈ [$kmin, $kmax]")
end

"""
    SyntheticFields.save_coefficients(fname::AbstractString, field::FourierBandVectorField)
    SyntheticFields.save_coefficients(g::Union{HDF5.File, HDF5.Group}, field::FourierBandVectorField)

Save Fourier coefficients to HDF5 file.

This can be useful for restarting a simulation keeping the same synthetic field.

See also [`SyntheticFields.load_coefficients!`](@ref).
"""
function save_coefficients(fname::AbstractString, field::FourierBandVectorField)
    endswith(".h5")(fname) || @warn lazy"saving HDF5 file: expected .h5 extension (got $fname)."
    h5open(ff -> save_coefficients(ff, field), fname, "w")
end

function save_coefficients(g::Union{HDF5.File, HDF5.Group}, field::FourierBandVectorField)
    (; qs, cs, Δks,) = field
    g["wavevectors_normalised"] = reinterpret(reshape, eltype(eltype(qs)), qs)
    g["wavevectors_step"] = collect(Δks)
    g["fourier_coefficients"] = reinterpret(reshape, eltype(eltype(cs)), cs)
    g
end

"""
    SyntheticFields.load_coefficients!(fname::AbstractString, field::FourierBandVectorField)
    SyntheticFields.load_coefficients!(g::Union{HDF5.File, HDF5.Group}, field::FourierBandVectorField)

Load Fourier coefficients from HDF5 file.

This can be useful for restarting a simulation keeping the same synthetic field.

See also [`SyntheticFields.save_coefficients`](@ref).
"""
function load_coefficients!(fname::AbstractString, field::FourierBandVectorField)
    endswith(".h5")(fname) || @warn lazy"loading HDF5 file: expected .h5 extension (got $fname)."
    h5open(ff -> load_coefficients!(ff, field), fname, "r")
end

function load_coefficients!(g::Union{HDF5.File, HDF5.Group}, field::FourierBandVectorField{T, N}) where {T, N}
    (; qs, cs, Δks,) = field
    # Check that Δks is the same as in the file. If that's not the case it's because the
    # domain sizes are not the same.
    Δks_file = read(g["wavevectors_step"]) :: Vector{T}
    length(Δks_file) == N || throw(DimensionMismatch("wrong dimensions of input file"))
    NTuple{N, T}(Δks_file) == Δks || error("domain periods from file are not the same as current ones")
    qs_file = read(g["wavevectors_normalised"]) :: Matrix{Int}
    cs_file = read(g["fourier_coefficients"])   :: Matrix{Complex{T}}
    @assert size(qs_file, 1) == size(cs_file, 1) == N
    Nc = size(qs_file, 2)
    @assert Nc == size(cs_file, 2)
    resize!(qs, Nc)
    resize!(cs, Nc)
    reinterpret(reshape, eltype(eltype(qs)), qs) .= qs_file
    reinterpret(reshape, eltype(eltype(cs)), cs) .= cs_file
    g
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

@kernel function from_fourier_grid_kernel!(
        op::F, cs::AbstractVector,
        @Const(qs::AbstractVector), @Const(Δks::NTuple{N}), @Const(ûs::NTuple{N})
    ) where {F, N}
    i = @index(Global, Linear)
    Ms = size(ûs[1])
    q⃗ = @inbounds qs[i]  # normalised wavevector
    js = ntuple(Val(N)) do d
        @inline
        q = q⃗[d]
        M = Ms[d]
        ifelse(q ≥ 0, q + 1, M + q + 1)  # source index
    end
    if all(js .< Ms)  # bounds check just in case (assumes one-based indexing)
        û = SVector(map(u -> @inbounds(u[js...]), ûs))
        k⃗ = SVector(q⃗ .* Δks)
        @inbounds cs[i] = @inline op(cs[i], û, k⃗)
    end
    nothing
end

"""
    SyntheticFields.from_fourier_grid!([op::Function], field::FourierBandVectorField{T, N}, ûs::NTuple{N, AbstractArray})

Copy values from vector field `ûs` in Fourier space onto a [`FourierBandVectorField`](@ref).

Only the wavenumbers within the band `[kmin, kmax]` in which `field` is defined are copied.

By default, old values in `field` are discarded and are simply replaced by the values in `ûs`.
One can use the optional `op` argument to change this behaviour, which should be a function
`op(old, new, k⃗)` taking 3 `SVector{N}`.
For example, passing `op(old, new, k⃗) = old - new * sum(abs2, k⃗)` means that new values are first multiplied by ``|\\bm{k}|^2``
and then subtracted from previously existent ones. The result is then written onto `field`.
"""
function from_fourier_grid!(op::F, field::FourierBandVectorField{T, N}, ûs::NTuple{N, AbstractArray}) where {F, T, N}
    (; qs, cs, Δks) = field
    Base.require_one_based_indexing(ûs...)  # assumed in GPU kernel
    backend = KA.get_backend(qs)
    groupsize = min(nextpow(2, length(qs)), 1024)
    ndrange = size(qs)
    kernel! = from_fourier_grid_kernel!(backend, groupsize, ndrange)
    kernel!(op, cs, qs, Δks, ûs)
    field
end

from_fourier_grid_default_op(old, new, k⃗) = new  # replace old values with new ones from Fourier-space fields
from_fourier_grid!(field::FourierBandVectorField, ûs::NTuple) = from_fourier_grid!(from_fourier_grid_default_op, field, ûs)

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
