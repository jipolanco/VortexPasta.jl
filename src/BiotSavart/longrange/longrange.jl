"""
    LongRangeCacheState

Describes the current state of the long-range fields in Fourier space.

It has two fields, which allow to know which field is currently stored and its smoothing scale
(width ``σ`` of Gaussian filter):

- `quantity::Symbol` (can be `:undef`, `:vorticity`, `:velocity`, `:streamfunction`)

- `smoothing_scale::Float64`: width ``σ`` of applied Gaussian filter. This is typically
  ``σ = 1 / α√2`` where ``α`` is Ewald's splitting parameter. Its value is 0 if the field
  has not been (yet) smoothed. This is usually the case after calling
  [`compute_vorticity_fourier!`](@ref).

See [`BiotSavart.get_longrange_field_fourier`](@ref) for more details.
"""
mutable struct LongRangeCacheState
    quantity :: Symbol  # quantity currently held by the cache (:undef, :vorticity, :velocity, :streamfunction)
    smoothing_scale :: Float64  # width σ of Gaussian filter (0 means unsmoothed)
end

LongRangeCacheState() = LongRangeCacheState(:undef, 0)

Base.copy(state::LongRangeCacheState) = LongRangeCacheState(state.quantity, state.smoothing_scale)

function Base.show(io::IO, state::LongRangeCacheState)
    (; quantity, smoothing_scale,) = state
    print(io, "LongRangeCacheState(quantity = $quantity, smoothing_scale = $smoothing_scale)")
    nothing
end

struct LongRangeCacheCommon{
        T <: AbstractFloat,
        Params <: ParamsLongRange,
        WaveNumbers <: NTuple{3, AbstractVector},
        PointCharges <: PointData{T},
        FourierVectorField <: StructArray{Vec3{Complex{T}}, 3},
        RealScalarField <: AbstractArray{T, 3},
        Timer <: TimerOutput,
    }
    # NOTE: the _d suffix means that data is on the device (i.e. GPU for GPU-based backends)
    params        :: Params
    wavenumbers_d :: WaveNumbers
    pointdata_d   :: PointCharges        # non-uniform data in physical space
    uhat_d        :: FourierVectorField  # uniform Fourier-space data (3 × [Nx, Ny, Nz])
    ewald_prefactor :: T                 # prefactor Γ/V (also included in ewald_op_d)
    ewald_op_d      :: RealScalarField    # real-valued Ewald operator in Fourier space ([Nx, Ny, Nz])
    state           :: LongRangeCacheState
    to_d            :: Timer             # timer for operations run on the device
end

function LongRangeCacheCommon(
        pcommon::ParamsCommon,
        params::ParamsLongRange,
        wavenumbers::NTuple{3, AbstractVector},
        pointdata::PointData{T},
        timer::TimerOutput,
    ) where {T}
    (; Γ, Ls, α,) = pcommon
    (; backend,) = params
    @assert T === eltype(pcommon)
    @assert α !== Zero()
    Nks = map(length, wavenumbers)
    uhat_d = init_fourier_vector_field(backend, T, Nks) :: StructArray{Vec3{Complex{T}}, 3}
    ewald_prefactor = Γ / prod(Ls)
    device = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    wavenumbers_d = adapt(device, wavenumbers)  # copy wavenumbers onto device if needed
    pointdata_d = adapt(device, pointdata)      # create PointData replica on the device if needed
    ewald_op_d = init_ewald_fourier_operator(T, backend, wavenumbers_d, α, ewald_prefactor)
    state = LongRangeCacheState()
    LongRangeCacheCommon(params, wavenumbers_d, pointdata_d, uhat_d, ewald_prefactor, ewald_op_d, state, timer)
end

has_real_to_complex(c::LongRangeCacheCommon) = has_real_to_complex(c.params)

# Initialise Fourier vector field with the right memory layout.
# This default implementation can be overridden by other backends.
# That's the case of FINUFFT, which needs a specific memory layout to perform simultaneous
# NUFFTs of the 3 vector field components.
function init_fourier_vector_field(backend::LongRangeBackend, ::Type{T}, Nks::Dims) where {T <: Real}
    device = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    # Initialising arrays in parallel (using KA) may be good for performance;
    # see https://juliagpu.github.io/KernelAbstractions.jl/dev/examples/numa_aware/.
    components = ntuple(Val(3)) do _
        KA.zeros(device, Complex{T}, Nks)
    end
    StructArray{Vec3{Complex{T}}}(components)
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
    (; common, backend,)= p
    # If long-range stuff is run on a GPU, we create a separate TimerOutput.
    # This is to avoid short- and long-range parts writing asynchronously to the same timer.
    to_d = KA.get_backend(backend) isa KA.CPU ? to : TimerOutput("Long range (GPU)")
    if common.α === Zero()
        NullLongRangeCache()  # disables Ewald method / long-range computations
    else
        init_cache_long_ewald(common, p, pointdata, to_d)
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

For convenience, point data already in `cache.common.pointdata_d` is copied to the new cache.
This means that, if one already filled the original cache using
[`add_point_charges!`](@ref), then one can directly call [`transform_to_fourier!`](@ref)
with the new cache to get the vorticity field in Fourier space at the wanted resolution `Ns`.
"""
function Base.similar(cache::LongRangeCache, Ns::Dims{3})
    params_old = cache.common.params :: ParamsLongRange
    (; backend, quad, common,) = params_old
    params_new = ParamsLongRange(backend, quad, common, Ns, params_old.truncate_spherical)
    pointdata = copy(cache.common.pointdata_d)
    BiotSavart.init_cache_long(params_new, pointdata)
end

"""
    transform_to_fourier!(cache::LongRangeCache)

Transform stored non-uniform data to Fourier space.

This usually corresponds to a type-1 NUFFT.

Non-uniform data must be first added via [`add_point_charges!`](@ref).
"""
function transform_to_fourier! end

@kernel function to_smoothed_streamfunction_kernel!(us::NTuple, @Const(ewald_op))
    I = @index(Global, Linear)
    for u ∈ us
        @inbounds u[I] *= ewald_op[I]
    end
    nothing
end

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
    (; uhat_d, state, ewald_op_d,) = c.common
    @assert size(uhat_d) === size(ewald_op_d)
    from_vorticity = state.quantity === :vorticity && state.smoothing_scale == 0
    @assert from_vorticity
    inds = eachindex(ewald_op_d, uhat_d)
    @assert inds isa AbstractUnitRange  # make sure we're using linear indexing (more efficient)
    ka_backend = KA.get_backend(c)
    kernel = ka_generate_kernel(to_smoothed_streamfunction_kernel!, ka_backend, uhat_d)
    uhat_comps = StructArrays.components(uhat_d)
    kernel(uhat_comps, ewald_op_d)
    state.quantity = :streamfunction
    state.smoothing_scale = ewald_smoothing_scale(c)
    uhat_d
end

# Applies Biot-Savart + Gaussian smoothing operators in Fourier space.
@kernel function velocity_from_vorticity_kernel!(uhat::NTuple, @Const(ks), @Const(ewald_op))
    I = @index(Global, Cartesian)
    @inbounds op = ewald_op[I]
    op_times_k⃗ = Vec3(map((k, i) -> @inbounds(op * k[i]), ks, Tuple(I)))
    u⃗ = Vec3(map(u -> @inbounds(u[I]), uhat))
    u⃗ = @inbounds im * u⃗
    # We use @noinline to avoid possible wrong results on the CPU!! (due to overoptimisation?)
    u⃗ = @noinline op_times_k⃗ × u⃗
    for n ∈ eachindex(uhat)
        @inbounds uhat[n][I] = u⃗[n]
    end
    nothing
end

# Applies curl operator in Fourier space.
@kernel function velocity_from_streamfunction_kernel!(uhat::NTuple, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), ks, Tuple(I)))
    u⃗ = Vec3(map(u -> @inbounds(u[I]), uhat))
    u⃗ = @inbounds im * u⃗
    # We use @noinline to avoid possible wrong results on the CPU!! (due to overoptimisation?)
    u⃗ = @noinline k⃗ × u⃗
    for n ∈ eachindex(uhat)
        @inbounds uhat[n][I] = u⃗[n]
    end
    nothing
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
    (; uhat_d, state, ewald_op_d, wavenumbers_d,) = c.common
    σ = ewald_smoothing_scale(c)
    @assert size(uhat_d) === size(ewald_op_d)
    from_vorticity = state.quantity === :vorticity && state.smoothing_scale == 0
    from_streamfunction = state.quantity === :streamfunction && state.smoothing_scale == σ
    @assert from_vorticity || from_streamfunction
    ka_backend = KA.get_backend(c)
    uhat_comps = StructArrays.components(uhat_d)
    if from_vorticity
        let kernel = ka_generate_kernel(velocity_from_vorticity_kernel!, ka_backend, uhat_d)
            kernel(uhat_comps, wavenumbers_d, ewald_op_d)
        end
    elseif from_streamfunction
        let kernel = ka_generate_kernel(velocity_from_streamfunction_kernel!, ka_backend, uhat_d)
            kernel(uhat_comps, wavenumbers_d)
        end
    end
    state.quantity = :velocity
    state.smoothing_scale = ewald_smoothing_scale(c)
    c
end

# These can be useful for generic code dealing with velocity, streamfunction or both.
to_smoothed_field!(::Velocity, c::LongRangeCache) = to_smoothed_velocity!(c)
to_smoothed_field!(::Streamfunction, c::LongRangeCache) = to_smoothed_streamfunction!(c)

"""
    interpolate_to_physical!([output::StructVector{<:Vec3},] cache::LongRangeCache) -> output

Perform type-2 NUFFT to interpolate values in `cache.common.uhat_d` to non-uniform
points in physical space.

Results are written to the `output` vector, which defaults to `cache.pointdata_d.charges`.
This vector is returned by this function.
"""
function interpolate_to_physical!(cache::LongRangeCache)
    output = cache.common.pointdata_d.charges
    interpolate_to_physical!(output, cache)
end

function interpolate_to_physical!(output, cache::LongRangeCache)
    # TODO: what if the output is in a different device? (CPU/GPU)?
    @assert typeof(output) === typeof(cache.common.pointdata_d.charges)
    _interpolate_to_physical!(output, cache)
    output
end

@kernel function init_ewald_fourier_operator_kernel!(
        u::AbstractArray{T, 3} where {T},
        @Const(ks),
        @Const(β::Real),  # = -1/(4α²)
        @Const(prefactor::Real),
    )
    I = @index(Global, Cartesian)
    k⃗ = map(getindex, ks, Tuple(I))
    k² = sum(abs2, k⃗)
    # Operator converting vorticity to coarse-grained streamfunction
    y = prefactor * exp(β * k²) / k²
    @inbounds u[I] = ifelse(iszero(k²), zero(y), y)
    nothing
end

function init_ewald_fourier_operator(
        ::Type{T}, backend::LongRangeBackend, ks, α::Real, prefactor::Real,
    ) where {T <: Real}
    ka_backend = KA.get_backend(backend)
    dims = map(length, ks)
    u = KA.zeros(ka_backend, T, dims)
    @assert adapt(ka_backend, ks) === ks  # check that `ks` vectors are already on the device
    β = -1 / (4 * α^2)
    kernel = ka_generate_kernel(init_ewald_fourier_operator_kernel!, ka_backend, u)
    kernel(u, ks, β, prefactor)
    u
end

function rescale_coordinates!(c::LongRangeCache)
    L = expected_period(backend(c))
    _rescale_coordinates!(c, L)
end

# Backend doesn't define `expected_period`, so no rescaling is needed.
_rescale_coordinates!(::LongRangeCache, ::Nothing) = nothing

function _rescale_coordinates!(c::LongRangeCache, L_expected::Real)
    (; Ls,) = c.common.params.common
    (; points,) = c.common.pointdata_d
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

@kernel function fold_coordinates_kernel!(points::NTuple, @Const(lims), @Const(L))
    i = @index(Global, Linear)
    for x ∈ points
        @inbounds x[i] = _fold_coordinate(x[i], lims, L)
    end
    nothing
end

function _fold_coordinates!(c::LongRangeCache, lims_in::NTuple{2, Real}, L_in::Real)
    (; points,) = c.common.pointdata_d
    points_comp = StructArrays.components(points) :: NTuple
    T = eltype(points_comp[1])
    @assert T <: AbstractFloat
    lims = convert.(T, lims_in)
    L = convert(T, L_in)
    ka_backend = KA.get_backend(c)
    @assert length(points) > 0
    # Instead of using ka_generate_kernel, we create a kernel with dynamically sized
    # ndrange, which means (I think) that the kernel won't be recompiled when length(points)
    # changes (which happens all the time during a simulation).
    workgroupsize = 1024
    kernel = fold_coordinates_kernel!(ka_backend, workgroupsize)
    kernel(points_comp, lims, L; ndrange = size(points))
    nothing
end

# We assume that L = b - a.
@inline function _fold_coordinate(x::T, (a, b)::NTuple{2, T}, L::T) where {T <: AbstractFloat}
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

The vorticity, written to `cache.common.uhat_d`, is estimated using some variant of
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
    (; uhat_d, state,) = cache.common
    rescale_coordinates!(cache)  # may be needed by the backend (e.g. FINUFFT requires period L = 2π)
    fold_coordinates!(cache)     # may be needed by the backend (e.g. FINUFFT requires x ∈ [-3π, 3π])
    transform_to_fourier!(cache)
    truncate_spherical!(cache)   # doesn't do anything if truncate_spherical = false (default)
    state.quantity = :vorticity
    state.smoothing_scale = 0
    uhat_d
end

@kernel function truncate_spherical_kernel!(uhat::NTuple, @Const(ks), @Const(kmax²))
    I = @index(Global, Cartesian)
    T = eltype(uhat[1])
    k⃗ = map((v, i) -> @inbounds(v[i]), ks, Tuple(I))
    k² = sum(abs2, k⃗)
    truncate = k² > kmax²
    for u ∈ uhat
        @inbounds u[I] = ifelse(truncate, zero(T), u[I])
    end
    nothing
end

function truncate_spherical!(cache)
    (; uhat_d, params, wavenumbers_d,) = cache.common
    if !params.truncate_spherical
        return  # don't apply spherical truncation
    end
    kmax = maximum_wavenumber(params)
    kmax² = kmax^2
    ka_backend = KA.get_backend(cache)
    uhat_comps = StructArrays.components(uhat_d)
    kernel = ka_generate_kernel(truncate_spherical_kernel!, ka_backend, uhat_d)
    kernel(uhat_comps, wavenumbers_d, kmax²)
    nothing
end

# Here pointdata_h is the point data on the CPU.
# It is used as an intermediate buffer in the GPU implementation.
# To avoid allocations, it should be passed by the caller. Its default value is mostly for
# backwards compatibility, when this function only accepted 2 arguments.
function set_interpolation_points!(
        cache::LongRangeCache,
        fs::VectorOfFilaments,
    )
    # Note that we only modify the `points` field of PointData. The other ones are not
    # needed for interpolation.
    (; pointdata_d,) = cache.common
    (; points, charges, points_h,) = pointdata_d
    Npoints = sum(length, fs)
    resize_no_copy!(points, Npoints)
    resize_no_copy!(charges, Npoints)  # this is the interpolation output
    ka_backend = KA.get_backend(cache)
    set_interpolation_points_impl!(ka_backend, fs, points, points_h)
    rescale_coordinates!(cache)
    fold_coordinates!(cache)
    nothing
end

# CPU implementation
function set_interpolation_points_impl!(
        ::KA.CPU, fs, points_d, points_h,
    )
    @assert isempty(points_h)  # never used in the CPU implementation (so it should stay empty)
    xs_d = StructArrays.components(points_d)  # (xs, ys, zs)
    _set_interpolation_points!(xs_d, fs)
    nothing
end

# GPU implementation: first write to pointdata_h (on the CPU), then copy to pointdata_d (on the GPU).
function set_interpolation_points_impl!(
        ka_backend::KA.GPU, fs, points_d, points_h,
    )
    Npoints = length(points_d)  # for now, only points_d has the right size
    resize_no_copy!(points_h, Npoints)  # resize temporary array (on the CPU)
    xs_d = StructArrays.components(points_d)  # (xs, ys, zs)
    xs_h = StructArrays.components(points_h)  # (xs, ys, zs)
    _set_interpolation_points!(xs_h, fs)  # gather points on the CPU
    # Copy to GPU
    for i ∈ eachindex(xs_d, xs_h)
        # KA.copyto!(ka_backend, xs_d[i], xs_h[i])  # may fail on CUDA due to pinning of CPU memory (https://github.com/JuliaGPU/CUDA.jl/issues/2594)
        copyto!(xs_d[i], xs_h[i])  # this doesn't fail (avoids pinning; probably slower but always works)
    end
    nothing
end

function _set_interpolation_points!(points::NTuple, fs)
    @assert KA.get_backend(points[1]) isa KA.CPU  # output is on the CPU
    Npoints = length(points[1])
    n = 0
    for f ∈ fs, X ∈ f
        n += 1
        for i ∈ eachindex(X)
            @inbounds points[i][n] = X[i]
        end
    end
    @assert n == Npoints
    points
end

# Here pointdata_h is used as a buffer in the GPU implementation.
# See set_interpolation_points! for details.
function add_long_range_output!(
        vs::AbstractVector{<:VectorOfVelocities}, cache::LongRangeCache,
        charges = cache.common.pointdata_d.charges,
    )
    (; charges_h,) = cache.common.pointdata_d
    nout = sum(length, vs)
    nout == length(charges) || throw(DimensionMismatch("wrong length of output vector `vs`"))
    ka_backend = KA.get_backend(cache)
    add_long_range_output_impl!(ka_backend, vs, charges, charges_h)
    vs
end

function add_long_range_output_impl!(::KA.CPU, vs, charges_d, charges_h)
    @assert isempty(charges_h)  # never used in the CPU implementation (so it should stay empty)
    _add_long_range_output!(vs, charges_d)
    nothing
end

# GPU implementation: first copy from charges_d (GPU) to charges_h (CPU), then add results
# to `vs`.
function add_long_range_output_impl!(ka_backend::KA.GPU, vs, charges_d, charges_h)
    # Device-to-host copy
    qs_d = StructArrays.components(charges_d)
    qs_h = StructArrays.components(charges_h)
    for i ∈ eachindex(qs_d, qs_h)
        resize_no_copy!(qs_h[i], length(qs_d[i]))
        # KA.copyto!(ka_backend, qs_h[i], qs_d[i])  # may fail on CUDA due to pinning of CPU memory (https://github.com/JuliaGPU/CUDA.jl/issues/2594)
        copyto!(qs_h[i], qs_d[i])  # this doesn't fail (avoids pinning; probably slower but always works)
    end
    KA.synchronize(ka_backend)  # make sure we're done copying data to CPU (needed on CUDA, where KA.copyto! is asynchronous)
    # Now add long-range values to `vs` output.
    _add_long_range_output!(vs, charges_h)
    nothing
end

function _add_long_range_output!(vs, charges::StructVector)
    n = 0
    @inbounds for v ∈ vs, j ∈ eachindex(v)
        q = charges[n += 1]
        v[j] = v[j] + real(q)
    end
    vs
end

function _ensure_hermitian_symmetry!(
        wavenumbers::NTuple{N},
        us::AbstractArray{<:Any, N},  # this can be a StructVector
    ) where {N}
    # Ensure Hermitian symmetry one dimension at a time.
    _ensure_hermitian_symmetry!(wavenumbers, Val(N), us)
end

@inline function _ensure_hermitian_symmetry!(wavenumbers, ::Val{d}, us) where {d}
    N = size(us, d)
    kd = wavenumbers[d]
    Δk = kd[2]
    is_r2c = last(kd) > 0  # real-to-complex transform
    T = eltype(us)
    if is_r2c
        @assert d == 1  # r2c transform is necessarily in the first dimension
        inds = ntuple(j -> j == 1 ? lastindex(kd) : Colon(), Val(ndims(us)))
        @inbounds @views us[inds...] .= Ref(zero(T))  # use Ref to avoid broadcasting over zero(T) (in case T <: Vec3)
    elseif iseven(N)  # case of complex-to-complex transform (nothing to do if N is odd)
        imin = (N ÷ 2) + 1  # asymmetric mode
        @assert -kd[imin] ≈ kd[imin - 1] + Δk
        inds = ntuple(j -> j == d ? imin : Colon(), Val(ndims(us)))
        @inbounds @views us[inds...] .= Ref(zero(T))
    end
    _ensure_hermitian_symmetry!(wavenumbers, Val(d - 1), us)
end

_ensure_hermitian_symmetry!(wavenumbers, ::Val{0}, us) = us  # we're done, do nothing

include("ka_utils.jl")
include("exact_sum.jl")
include("finufft.jl")
include("nonuniformffts.jl")
