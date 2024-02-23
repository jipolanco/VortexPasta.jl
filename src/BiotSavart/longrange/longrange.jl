mutable struct LongRangeCacheState
    quantity :: Symbol  # quantity currently held by the cache (:undef, :vorticity, :velocity, :streamfunction)
    smoothed :: Bool    # true if Ewald's Gaussian filter has already been applied
end

LongRangeCacheState() = LongRangeCacheState(:undef, false)

struct LongRangeCacheCommon{
        T <: AbstractFloat,
        Params <: ParamsLongRange,
        WaveNumbers <: NTuple{3, AbstractVector},
        PointCharges <: PointData{T},
        FourierVectorField <: StructArray{Vec3{Complex{T}}, 3},
        Timer <: TimerOutput,
    }
    params      :: Params
    wavenumbers :: WaveNumbers
    pointdata   :: PointCharges        # non-uniform data in physical space
    uhat        :: FourierVectorField  # uniform Fourier-space data (3 × [Nx, Ny, Nz])
    ewald_prefactor :: T               # prefactor Γ/V (also included in ewald_op)
    ewald_op    :: Array{T, 3}         # Ewald operator in Fourier space ([Nx, Ny, Nz])
    state       :: LongRangeCacheState
    to          :: Timer
end

function LongRangeCacheCommon(
        pcommon::ParamsCommon,
        params::ParamsLongRange,
        wavenumbers::NTuple{3, AbstractVector},
        pointdata::PointData{T},
        timer::TimerOutput,
    ) where {T}
    (; Γ, Ls, α,) = pcommon
    @assert T === eltype(pcommon)
    @assert α !== Zero()
    Nks = map(length, wavenumbers)
    uhat = StructArray{Vec3{Complex{T}}}(undef, Nks)
    ewald_prefactor = Γ / prod(Ls)
    ewald_op = init_ewald_fourier_operator(T, wavenumbers, α, ewald_prefactor)
    state = LongRangeCacheState()
    LongRangeCacheCommon(params, wavenumbers, pointdata, uhat, ewald_prefactor, ewald_op, state, timer)
end

"""
    NullLongRangeCache <: LongRangeCache

Dummy cache type returned by [`init_cache_long`](@ref) when long-range
computations are disabled.

This is the case when the Ewald splitting parameter ``α`` is set to `Zero()`.
"""
struct NullLongRangeCache <: LongRangeCache end

backend(::NullLongRangeCache) = NullLongRangeBackend()

"""
    init_cache_long(p::ParamsLongRange, pointdata::PointData, [to::TimerOutput]) -> LongRangeCache

Initialise the cache for the long-range backend defined in `p`.

Note that, if `pc.α === Zero()`, then long-range computations are disabled and
this returns a [`NullLongRangeCache`](@ref).
"""
function init_cache_long(p::ParamsLongRange, pointdata::PointData, to = TimerOutput())
    pc = p.common
    if pc.α === Zero()
        NullLongRangeCache()  # disables Ewald method / long-range computations
    else
        init_cache_long_ewald(pc, p, pointdata, to)
    end
end

"""
    Base.similar(cache::LongRangeCache, Ns::Dims{3}) -> LongRangeCache

Create a [`LongRangeCache`](@ref) similar to an existent one, but possibly holding a
different amount of wavenumbers in Fourier space.

In principle, the grid resolution `Ns` can be different from the original one
(`cache.common.params.Ns`).
This can be useful for computing high-resolution fields in Fourier space (e.g. for extending
the data to higher wavenumbers than allowed by the original cache).

For convenience, point data already in `cache.common.pointdata` is copied to the new cache.
This means that, if one already filled the original cache using
[`add_point_charges!`](@ref), then one can directly call [`transform_to_fourier!`](@ref)
with the new cache to get the vorticity field in Fourier space at the wanted resolution `Ns`.
"""
function Base.similar(cache::LongRangeCache, Ns::Dims{3})
    params_old = cache.common.params :: ParamsLongRange
    (; backend, quad, common,) = params_old
    params_new = ParamsLongRange(backend, quad, common, Ns)
    pointdata = copy(cache.common.pointdata)
    BiotSavart.init_cache_long(params_new, pointdata)
end

"""
    transform_to_fourier!(cache::LongRangeCache)

Transform stored non-uniform data to Fourier space.

This usually corresponds to a type-1 NUFFT.

Non-uniform data must be first added via [`add_point_charges!`](@ref).
"""
function transform_to_fourier! end

@doc raw"""
    to_smoothed_streamfunction!(cache::LongRangeCache)

Convert Fourier-transformed vorticity field to coarse-grained streamfunction field in
Fourier space.

This operation can be simply written as:

```math
\hat{\bm{ψ}}_{α}(\bm{k}) = ϕ_α(|\bm{k}|) \, \hat{\bm{ω}}(\bm{k})
```

where

```math
ϕ_α(k) = \frac{Γ}{V} \, \frac{e^{-k^2 / 4α^2}}{k^2}
```

is the Ewald operator. The effect of this operator is to:

1. invert the Laplacian in ``-∇² \bm{ψ} = \bm{\omega}``;
2. smoothen the fields according to the Ewald parameter ``α`` (an inverse length scale);
3. rescale values by the vortex circulation ``Γ`` and the volume ``V`` of a periodic cell
   so that the streamfunction has the right units (``L^2 T^{-1}``).

This function should be called after [`compute_vorticity_fourier!`](@ref).
If one only needs the velocity and not the streamfunction, one can also directly call
[`to_smoothed_velocity!`](@ref).
"""
function to_smoothed_streamfunction!(c::LongRangeCache)
    (; uhat, state, ewald_op,) = c.common
    @assert size(uhat) === size(ewald_op)
    from_vorticity = state.quantity === :vorticity && !state.smoothed
    @assert from_vorticity
    inds = eachindex(ewald_op, uhat)
    @assert inds isa AbstractUnitRange  # make sure we're using linear indexing (more efficient)
    @inbounds for i ∈ inds
        uhat[i] = ewald_op[i] * uhat[i]
    end
    state.quantity = :streamfunction
    state.smoothed = true
    uhat
end

@doc raw"""
    to_smoothed_velocity!(cache::LongRangeCache)

Convert Fourier-transformed vorticity field to coarse-grained velocity field in
Fourier space.

If called right after [`compute_vorticity_fourier!`](@ref), this function performs the
operation:

```math
\hat{\bm{v}}_{α}(\bm{k}) = i \bm{k} × ϕ_α(|\bm{k}|) \, \, \hat{\bm{ω}}(\bm{k}),
```

where ``ϕ_α`` is the Ewald operator defined in [`to_smoothed_streamfunction!`](@ref).

Optionally, if one is also interested in the streamfunction, one can call
[`to_smoothed_streamfunction!`](@ref) *before* this function.
In that case, the cache already contains the smoothed streamfunction, and only the curl
operator (``i \bm{k} ×``) is applied by this function.
"""
function to_smoothed_velocity!(c::LongRangeCache)
    (; uhat, state, ewald_op, wavenumbers,) = c.common
    @assert size(uhat) === size(ewald_op)
    from_vorticity = state.quantity === :vorticity && !state.smoothed
    from_streamfunction = state.quantity === :streamfunction && state.smoothed
    @assert from_vorticity || from_streamfunction
    if from_vorticity
        @inbounds for I ∈ CartesianIndices(ewald_op)
            op = ewald_op[I]
            op_times_k⃗ = Vec3(map((k, i) -> @inbounds(op * k[i]), wavenumbers, Tuple(I)))
            uhat[I] = op_times_k⃗ × (im * uhat[I])
        end
    elseif from_streamfunction
        @inbounds for I ∈ CartesianIndices(ewald_op)
            k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), wavenumbers, Tuple(I)))
            uhat[I] = k⃗ × (im * uhat[I])
        end
    end
    state.quantity = :velocity
    state.smoothed = true
    c
end

"""
    interpolate_to_physical!(cache::LongRangeCache)

Perform type-2 NUFFT to interpolate values in `cache.uhat` to non-uniform
points in physical space.

Results are written to `cache.pointdata.charges`.
"""
function interpolate_to_physical! end

function init_ewald_fourier_operator!(
        u::AbstractArray{T, 3} where {T}, ks, α::Real, prefactor::Real,
    )
    β = -1 / (4 * α^2)
    for I ∈ CartesianIndices(u)
        k⃗ = map(getindex, ks, Tuple(I))
        k² = sum(abs2, k⃗)
        # Operator converting vorticity to coarse-grained streamfunction
        y = prefactor * exp(β * k²) / k²
        u[I] = ifelse(iszero(k²), zero(y), y)
    end
    u
end

function init_ewald_fourier_operator(::Type{T}, ks, args...) where {T <: Real}
    dims = map(length, ks)
    u = Array{T}(undef, dims)
    init_ewald_fourier_operator!(u, ks, args...)
end

function rescale_coordinates!(c::LongRangeCache)
    L = expected_period(backend(c))
    _rescale_coordinates!(c, L)
end

# Backend doesn't define `expected_period`, so no rescaling is needed.
_rescale_coordinates!(::LongRangeCache, ::Nothing) = nothing

function _rescale_coordinates!(c::LongRangeCache, L_expected::Real)
    (; Ls,) = c.common.params.common
    (; points,) = c.common.pointdata
    for (xs, L) ∈ zip(StructArrays.components(points), Ls)
        _rescale_coordinates!(xs, L, L_expected)
    end
    nothing
end

function _rescale_coordinates!(xs::AbstractVector, L::Real, L_expected)
    if L != L_expected
        xs .*= L_expected / L
    end
    nothing
end

# Note: This function must be called **after** `rescale_coordinates!`.
function fold_coordinates!(c::LongRangeCache)
    lims = folding_limits(backend(c))
    L = expected_period(backend(c))
    _fold_coordinates!(c, lims, L)
end

# Backend doesn't define `folding_limits`, so no rescaling is needed.
_fold_coordinates!(::LongRangeCache, ::Nothing, ::Any) = nothing

function _fold_coordinates!(c::LongRangeCache, lims::NTuple{2}, L::Real)
    for xs ∈ StructArrays.components(c.common.pointdata.points)
        @inbounds for (i, x) ∈ pairs(xs)
            xs[i] = _fold_coordinate(x, lims, L)
        end
    end
    nothing
end

# We assume that L = b - a.
@inline function _fold_coordinate(x::Real, (a, b), L::Real)
    while x ≥ b
        x -= L
    end
    while x < a
        x += L
    end
    x
end

# Determines whether the backend requires non-uniform values to be complex.
# By default this is not the case, but backends needing that (such as FINUFFT) return
# Complex{T} instead of T.
non_uniform_type(::Type{T}, ::LongRangeBackend) where {T <: AbstractFloat} = T

"""
    compute_vorticity_fourier!(cache::LongRangeCache)

Estimate vorticity in Fourier space.

The vorticity, written to `cache.common.uhat`, is estimated using some variant of
non-uniform Fourier transforms (depending on the chosen backend).

Note that this function doesn't perform smoothing over the vorticity using the Ewald operator.
Moreover, the resulting vorticity doesn't have the right dimensions, as it must be
multiplied by ``Γ/V`` (where ``Γ`` is the circulation and ``V`` is the volume of a periodic
cell) to have dimensions ``T^{-1}``. In fact, this factor is included in the Ewald operator
(see [`to_smoothed_streamfunction!`](@ref) for details).

Must be called after [`add_point_charges!`](@ref).

After calling this function, one may want to use [`to_smoothed_streamfunction!`](@ref) and/or
[`to_smoothed_velocity!`](@ref) to obtain the respective fields.
"""
function compute_vorticity_fourier!(cache::LongRangeCache)
    (; uhat, state,) = cache.common
    rescale_coordinates!(cache)  # may be needed by the backend (e.g. FINUFFT requires period L = 2π)
    fold_coordinates!(cache)     # may be needed by the backend (e.g. FINUFFT requires x ∈ [-3π, 3π])
    transform_to_fourier!(cache)
    state.quantity = :vorticity
    state.smoothed = false
    uhat
end

function set_interpolation_points!(
        cache::LongRangeCache,
        fs::VectorOfFilaments,
    )
    (; pointdata,) = cache.common
    Npoints = sum(length, fs)
    set_num_points!(pointdata, Npoints)
    n = 0
    for f ∈ fs, X ∈ f
        add_point!(pointdata, X, n += 1)
    end
    @assert n == Npoints
    rescale_coordinates!(cache)
    fold_coordinates!(cache)
    nothing
end

function add_long_range_output!(
        vs::AbstractVector{<:VectorOfVelocities}, cache::LongRangeCache,
    )
    (; charges,) = cache.common.pointdata
    nout = sum(length, vs)
    nout == length(charges) || throw(DimensionMismatch("wrong length of output vector `vs`"))
    n = 0
    @inbounds for v ∈ vs, i ∈ eachindex(v)
        q = charges[n += 1]
        v[i] = v[i] + real(q)
    end
    vs
end

function _ensure_hermitian_symmetry!(wavenumbers::NTuple{N}, us::Array{<:Complex, N}) where {N}
    # Ensure Hermitian symmetry one dimension at a time.
    _ensure_hermitian_symmetry!(wavenumbers, Val(N), us)
end

@inline function _ensure_hermitian_symmetry!(wavenumbers, ::Val{d}, us) where {d}
    N = size(us, d)
    kd = wavenumbers[d]
    Δk = kd[2]
    is_r2c = last(kd) > 0  # real-to-complex transform
    if is_r2c
        @assert d == 1  # r2c transform is necessarily in the first dimension
        inds = ntuple(j -> j == 1 ? lastindex(kd) : Colon(), Val(ndims(us)))
        @inbounds @views us[inds...] .= 0
    elseif iseven(N)  # case of complex-to-complex transform (nothing to do if N is odd)
        imin = (N ÷ 2) + 1  # asymmetric mode
        @assert -kd[imin] ≈ kd[imin - 1] + Δk
        inds = ntuple(j -> j == d ? imin : Colon(), Val(ndims(us)))
        @inbounds @views us[inds...] .= 0
    end
    _ensure_hermitian_symmetry!(wavenumbers, Val(d - 1), us)
end

_ensure_hermitian_symmetry!(wavenumbers, ::Val{0}, us) = us  # we're done, do nothing

include("exact_sum.jl")
include("finufft.jl")
include("nonuniformffts.jl")
