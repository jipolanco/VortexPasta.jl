include("cache_common.jl")

# Initialise Fourier vector field with the right memory layout.
# This default implementation can be overridden by other backends.
# That's the case of the old FINUFFT backend (which has been removed), which needs a
# specific memory layout to perform simultaneous NUFFTs of the 3 vector field components.
function init_fourier_vector_field(backend::LongRangeBackend, ::Type{T}, Nks::Dims) where {T <: Real}
    ka_backend = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    # Initialising arrays in parallel (using KA) may be good for performance;
    # see https://juliagpu.github.io/KernelAbstractions.jl/dev/examples/numa_aware/.
    components = ntuple(Val(3)) do _
        KA.zeros(ka_backend, Complex{T}, Nks)
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
    init_cache_long(p::ParamsBiotSavart, pointdata::PointData) -> LongRangeCache

Initialise the cache for the long-range backend defined in `p.longrange`.

Note that, if `p.α === Zero()`, then long-range computations are disabled and
this returns a [`NullLongRangeCache`](@ref).
"""
function init_cache_long(params::ParamsBiotSavart, pointdata::PointData)
    # If long-range stuff is run on a GPU, we create a separate TimerOutput.
    # This is to avoid short- and long-range parts writing asynchronously to the same timer.
    if params.α === Zero()
        NullLongRangeCache()  # disables Ewald method / long-range computations
    else
        init_cache_long_ewald(params, params.longrange, pointdata)
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
    params_new = ParamsLongRange(backend, quad, common, Ns, params_old.truncate_spherical)
    params_bs = let p = get_parameters(cache)::ParamsBiotSavart
        ParamsBiotSavart(p.common, p.shortrange, params_new)
    end
    pointdata = copy(cache.common.pointdata)
    BiotSavart.init_cache_long(params_bs, pointdata)
end

"""
    transform_to_fourier!(cache::LongRangeCache, prefactor::Real)

Transform stored non-uniform data to Fourier space.

This usually corresponds to a type-1 NUFFT.

The `prefactor` is a real value that will be used to multiply each non-uniform value (or
equivalently each uniform value in Fourier space).

Non-uniform data must be first added via [`add_point_charges!`](@ref).
"""
function transform_to_fourier! end

@doc raw"""
    compute_streamfunction_fourier!(cache::LongRangeCache)

Convert Fourier-transformed vorticity field to streamfunction field in Fourier space.

This simply corresponds to inverting a Laplacian (``-∇^2 \bm{ψ} = \bm{ω}``), which in
Fourier space is:

```math
\hat{\bm{ψ}}_{α}(\bm{k}) = \hat{\bm{ω}}(\bm{k}) / k²
```

This function should be called after [`compute_vorticity_fourier!`](@ref).
If one only needs the velocity and not the streamfunction, one can also directly call
[`compute_velocity_fourier!`](@ref).
"""
function compute_streamfunction_fourier!(c::LongRangeCache)
    (; uhat, state,) = c.common
    wavenumbers = get_wavenumbers(c)
    from_vorticity = state.quantity === :vorticity && state.smoothing_scale == 0
    @assert from_vorticity
    ka_backend = KA.get_backend(c)
    kernel = ka_generate_kernel(fourier_inverse_laplacian_kernel!, ka_backend, uhat)
    uhat_comps = StructArrays.components(uhat)
    kernel(uhat_comps, wavenumbers)
    state.quantity = :streamfunction
    uhat
end

# This inverts -∇²ψ = ω  <-->  k² ψ̂ = ω̂
@kernel function fourier_inverse_laplacian_kernel!(us::NTuple, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), ks, Tuple(I)))
    k² = sum(abs2, k⃗)
    k²_inv = ifelse(iszero(k²), zero(k²), inv(k²))  # this implicitly sets the k⃗ = 0 mode to zero
    for u ∈ us
        @inbounds u[I] *= k²_inv
    end
    nothing
end

# Applies Biot-Savart + Gaussian smoothing operators in Fourier space.
@kernel function fourier_biot_savart_kernel!(uhat::NTuple, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), ks, Tuple(I)))
    u⃗ = Vec3(map(u -> @inbounds(u[I]), uhat))
    u⃗ = @inbounds im * u⃗
    # We use @noinline to avoid possible wrong results on the CPU!! (due to overoptimisation?) -- XXX: is this still the case?
    u⃗ = @noinline k⃗ × u⃗
    k² = sum(abs2, k⃗)
    k²_inv = ifelse(iszero(k²), zero(k²), inv(k²))  # this implicitly sets the k⃗ = 0 mode to zero
    u⃗ = u⃗ * k²_inv
    for n ∈ eachindex(uhat)
        @inbounds uhat[n][I] = u⃗[n]
    end
    nothing
end

# Applies curl operator in Fourier space.
@kernel function fourier_curl_kernel!(uhat::NTuple, @Const(ks))
    I = @index(Global, Cartesian)
    k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), ks, Tuple(I)))
    u⃗ = Vec3(map(u -> @inbounds(u[I]), uhat))
    u⃗ = @inbounds im * u⃗
    # We use @noinline to avoid possible wrong results on the CPU!! (due to overoptimisation?) -- XXX: is this still the case?
    u⃗ = @noinline k⃗ × u⃗
    for n ∈ eachindex(uhat)
        @inbounds uhat[n][I] = u⃗[n]
    end
    nothing
end

# Applies curl operator in Fourier space. The input is usually the velocity field smoothed
# via Ewald's operator (with smoothing scale ℓ_old = 1 / α√2).
@kernel function coarse_grained_curl_kernel!(uhat::NTuple, @Const(ks), @Const(ℓ::Real), @Const(ℓ_old::Real))
    I = @index(Global, Cartesian)
    k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), ks, Tuple(I)))
    k² = sum(abs2, k⃗)
    op = exp(-k² * (ℓ^2 - ℓ_old^2) / 2)
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

# This kernel just applies a Gaussian filter without any differential operators.
@kernel function coarse_graining_kernel!(uhat::NTuple, @Const(ks), @Const(ℓ::Real), @Const(ℓ_old::Real))
    I = @index(Global, Cartesian)
    k⃗ = Vec3(map((k, i) -> @inbounds(k[i]), ks, Tuple(I)))
    k² = sum(abs2, k⃗)
    op = exp(-k² * (ℓ^2 - ℓ_old^2) / 2)
    for n ∈ eachindex(uhat)
        @inbounds uhat[n][I] *= op
    end
    nothing
end

@doc raw"""
    compute_velocity_fourier!(cache::LongRangeCache)

Convert Fourier-transformed vorticity field to velocity field in Fourier space.

If called right after [`compute_vorticity_fourier!`](@ref), this function performs the
operation:

```math
\hat{\bm{v}}_{α}(\bm{k}) = \frac{i \bm{k} × \, \hat{\bm{ω}}(\bm{k})}{k^2}.
```

Optionally, if one is also interested in the streamfunction, one can call
[`compute_streamfunction_fourier!`](@ref) *before* this function.
In that case, the cache already contains the streamfunction, and only the curl
operator (``i \bm{k} ×``) is applied by this function.
"""
function compute_velocity_fourier!(c::LongRangeCache)
    (; uhat, state,) = c.common
    wavenumbers = get_wavenumbers(c)
    from_vorticity = state.quantity === :vorticity && state.smoothing_scale == 0
    from_streamfunction = state.quantity === :streamfunction && state.smoothing_scale == 0
    @assert from_vorticity || from_streamfunction
    ka_backend = KA.get_backend(c)
    uhat_comps = StructArrays.components(uhat)
    if from_vorticity
        let kernel = ka_generate_kernel(fourier_biot_savart_kernel!, ka_backend, uhat)
            kernel(uhat_comps, wavenumbers)
        end
    elseif from_streamfunction
        let kernel = ka_generate_kernel(fourier_curl_kernel!, ka_backend, uhat)
            kernel(uhat_comps, wavenumbers)
        end
    end
    state.quantity = :velocity
    c
end

"""
    to_coarse_grained_vorticity!(c::LongRangeCache, ℓ::Real) -> NamedTuple

Compute coarse-grained vorticity in Fourier space.

The vorticity is coarse-grained at scale ``ℓ`` using a Gaussian filter.

Ideally, for accuracy reasons, the smoothing scale ``ℓ`` should be larger than the Ewald
splitting scale ``σ = 1 / α√2``. Also note that applying this function leads to **loss of
small-scale information**, so that one should be careful when performing additional
calculations (energy spectra, ...) afterwards.

This can be useful for visualisations or physical analysis. It is _not_ used to compute
Biot–Savart velocities, and the scale ``ℓ`` need not be similar to the splitting lengthscale
in Ewald's method.

This function also sets the [`LongRangeCacheState`](@ref) to `quantity = :vorticity` and `smoothing_scale = ℓ`.

This function returns the same thing as [`get_longrange_field_fourier`](@ref), including the
resulting coarse-grained vorticity in Fourier space.

## Example

From a set of vortex filaments, evaluate their self-induced velocity (using
[`compute_on_nodes!`](@ref)) and then the resulting coarse-grained vorticity at their
locations:

```julia
# First compute the velocity on the filament nodes.
# As an intermediate result, the velocity in Fourier space is stored in the cache.
vs = map(similar ∘ nodes, fs)  # one velocity vector per filament node
fields = (; velocity = vs,)
cache = BiotSavart.init_cache(params)
compute_on_nodes!(fields, cache, fs)

# Now compute the coarse-grained vorticity (result is stored in the cache)
ℓ = 0.1  # smoothing scale
BiotSavart.to_coarse_grained_vorticity!(cache.longrange, ℓ)

# Finally, evaluate the resulting vorticity on filament nodes
ωs_ℓ = map(similar, vs)
BiotSavart.interpolate_to_physical!(cache.longrange)       # interpolate to vortex positions (result stored in the cache)
BiotSavart.copy_long_range_output!(ωs_ℓ, cache.longrange)  # copy results from cache to output array
```
"""
function to_coarse_grained_vorticity!(c::LongRangeCache, ℓ::Real; warn = true)
    (; uhat, state,) = c.common
    wavenumbers = get_wavenumbers(c)
    σ_ewald = ewald_smoothing_scale(c)
    if warn && ℓ < σ_ewald
        @warn lazy"for accuracy reasons, the smoothing scale (ℓ = $ℓ) should be larger than the Ewald splitting scale (σ = $σ_ewald)"
    end
    ℓ_old = state.smoothing_scale  # current smoothing scale (usually 0, but can be different)
    ka_backend = KA.get_backend(c)
    uhat_comps = StructArrays.components(uhat)
    if state.quantity === :velocity
        let kernel = ka_generate_kernel(coarse_grained_curl_kernel!, ka_backend, uhat)
            kernel(uhat_comps, wavenumbers, ℓ, ℓ_old)
        end
    elseif state.quantity === :vorticity
        let kernel = ka_generate_kernel(coarse_graining_kernel!, ka_backend, uhat)
            kernel(uhat_comps, wavenumbers, ℓ, ℓ_old)
        end
    else
        error(lazy"expected either the velocity or vorticity to be in the cache (got state.quantity = $(state.quantity))")
    end
    state.quantity = :vorticity
    state.smoothing_scale = ℓ
    get_longrange_field_fourier(c)
end

# These can be useful for generic code dealing with velocity, streamfunction or both.
compute_field_fourier!(::Velocity, c::LongRangeCache) = compute_velocity_fourier!(c)
compute_field_fourier!(::Streamfunction, c::LongRangeCache) = compute_streamfunction_fourier!(c)

"""
    interpolate_to_physical!([callback::Function], [output::StructVector{<:Vec3}], cache::LongRangeCache) -> output

Perform type-2 NUFFT to interpolate values in `cache.common.uhat` to non-uniform
points in physical space.

Results are written to the `output` vector, which defaults to `cache.outputs[1]`.
This vector is returned by this function.

## Using callbacks

The optional `callback` function should be a function `f(û, idx)` which takes:

- `û::NTuple{3, <:Complex}` a Fourier coefficient;
- `idx::NTuple{3, Int}` the index of `û` on the Fourier grid.

The function should return a tuple with the same type as `û`.

This can be used to modify the Fourier-space fields to be interpolated, but without really
modifying the values in `uhat` (and thus the state of the cache). For example, to
apply a Gaussian filter before interpolating:

```julia
σ = 0.1  # width of Gaussian filter (in physical space)
ks = BiotSavart.get_wavenumbers(cache.longrange)

callback(û::NTuple{3}, idx::NTuple{3}) = let
    k⃗ = @inbounds getindex.(ks, idx)
    k² = sum(abs2, k⃗)
    û .* exp(-σ^2 * k² / 2)
end

interpolate_to_physical!(callback, cache)
```

!!! note "GPU usage"

    As written above, the GPU kernels using the callback might not compile correctly.
    If that is the case, a possible solution is to use a callable struct instead.
    For example:

    ```julia
    struct GaussianFilter{WaveNumbers} <: Function
        α::Float64
        ks::WaveNumbers
    end

    @inline function (f::GaussianFilter)(û::NTuple{3}, idx::NTuple{3})
        (; α, ks) = f
        k⃗ = @inbounds getindex.(ks, idx)
        k² = sum(abs2, k⃗)
        û .* exp(-σ^2 * k² / 2)
    end

    callback = GaussianFilter(0.1, BiotSavart.get_wavenumbers(cache.longrange))  # this is now a callable object
    interpolate_to_physical(callback, cache)
    ```

    If this still doesn't work, one might need to replace broadcasts (`.`) with calls to `map` or `ntuple`.

"""
function interpolate_to_physical! end

@inline default_callback_interp(û::NTuple{3}, idx::NTuple{3}) = û  # by default we leave the original value unmodified

interpolate_to_physical!(cache::LongRangeCache) = interpolate_to_physical!(default_callback_interp, cache)
interpolate_to_physical!(output::StructVector, cache::LongRangeCache) = interpolate_to_physical!(default_callback_interp, output, cache)

function interpolate_to_physical!(callback::F, cache::LongRangeCache) where {F}
    output = default_interpolation_output(cache)
    resize_no_copy!(output, length(cache.pointdata.nodes))
    interpolate_to_physical!(callback, output, cache)
end

function interpolate_to_physical!(callback::F, output::StructVector, cache::LongRangeCache) where {F}
    # TODO: what if the output is in a different device? (CPU/GPU)?
    @assert typeof(output) === typeof(cache.common.pointdata.charges)
    _interpolate_to_physical!(callback, output, cache)
    output
end

@kernel function init_ewald_gaussian_operator_kernel!(
        u::AbstractArray{T, 3} where {T},
        @Const(ks),
        @Const(beta::Real),  # = -1/(4α²) // CUDA doesn't like Unicode characters?? ("ptxas fatal : Unexpected non-ASCII character encountered...")
    )
    I = @index(Global, Cartesian)
    k⃗ = map(getindex, ks, Tuple(I))
    k² = sum(abs2, k⃗)
    # This is simply a Gaussian smoothing operator in Fourier space (note: β = -1/4α²).
    @inbounds u[I] = exp(beta * k²)
    nothing
end

function init_ewald_gaussian_operator(
        ::Type{T}, backend::LongRangeBackend, ks, α::Real,
    ) where {T <: Real}
    ka_backend = KA.get_backend(backend)
    dims = map(length, ks)
    u = KA.zeros(ka_backend, T, dims)
    @assert adapt(ka_backend, ks) === ks  # check that `ks` vectors are already on the device
    β = -1 / (4 * α^2)
    kernel = ka_generate_kernel(init_ewald_gaussian_operator_kernel!, ka_backend, u)
    kernel(u, ks, β)
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
    (; points, nodes) = c.common.pointdata
    for (xs, L) ∈ zip(StructArrays.components(points), Ls)
        _rescale_coordinates!(xs, L, L_expected)
    end
    for (xs, L) ∈ zip(StructArrays.components(nodes), Ls)
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
        @inbounds x[i] = _fold_coordinate_between_lims(x[i], lims, L)
    end
    nothing
end

function _fold_coordinates!(c::LongRangeCache, lims_in::NTuple{2, Real}, L_in::Real)
    (; points, nodes) = c.common.pointdata
    points_comp = StructArrays.components(points) :: NTuple
    nodes_comp = StructArrays.components(nodes) :: NTuple
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
    kernel(nodes_comp, lims, L; ndrange = size(nodes))
    nothing
end

# We assume that L = b - a.
@inline function _fold_coordinate_between_lims(x::T, (a, b)::NTuple{2, T}, L::T) where {T <: AbstractFloat}
    while x ≥ b
        x -= L
    end
    while x < a
        x += L
    end
    x
end

"""
    process_point_charges!(cache::LongRangeCache)

Process list of point charges in long-range cache.

For long-range computations, this should be called after quadrature points and interpolation
nodes have been set, either via `add_point_charges!(cache, fs)` or by directly modifying
`cache.pointdata`. It must be called before any calls to
[`compute_vorticity_fourier!`](@ref) or [`interpolate_to_physical!`](@ref).

This function will process (and possibly modify) data in `cache.pointdata`. Therefore,
**it should only be called once** after points have been set. For example, the
[`NonuniformFFTsBackend`](@ref) assumes a domain of size `[0, 2π]ᵈ`, and thus a point
transformation is needed if the domain has a different size.
"""
function process_point_charges!(cache::LongRangeCache)
    rescale_coordinates!(cache)  # may be needed by the backend (e.g. NonuniformFFTs.jl requires period L = 2π)
    fold_coordinates!(cache)     # may be needed by the backend (e.g. NonuniformFFTs.jl prefers x ∈ [0, 2π])
    cache
end

"""
    compute_vorticity_fourier!(cache::LongRangeCache)

Estimate vorticity in Fourier space.

The vorticity, written to `cache.common.uhat`, is estimated using some kind of
non-uniform Fourier transforms (depending on the chosen backend).

Note that this function doesn't perform smoothing over the vorticity using the Ewald operator.

Must be called after [`add_point_charges!`](@ref).

After calling this function, one may want to use [`compute_streamfunction_fourier!`](@ref) and/or
[`compute_velocity_fourier!`](@ref) to obtain the respective fields.
"""
function compute_vorticity_fourier!(cache::LongRangeCache)
    (; uhat, state, vorticity_prefactor,) = cache.common
    transform_to_fourier!(cache, vorticity_prefactor)
    truncate_spherical!(cache)   # doesn't do anything if truncate_spherical = false (default)
    state.quantity = :vorticity
    state.smoothing_scale = 0
    uhat
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
    (; uhat, params,) = cache.common
    wavenumbers = get_wavenumbers(cache)
    if !params.truncate_spherical
        return  # don't apply spherical truncation
    end
    kmax = maximum_wavenumber(params)
    kmax² = kmax^2
    ka_backend = KA.get_backend(cache)
    uhat_comps = StructArrays.components(uhat)
    kernel = ka_generate_kernel(truncate_spherical_kernel!, ka_backend, uhat)
    kernel(uhat_comps, wavenumbers, kmax²)
    nothing
end

function set_interpolation_points!(cache::LongRangeCache, fs::VectorOfFilaments)
    @warn(
        """
        The BiotSavart.set_interpolation_points!(cache, fs) function is no longer needed and is deprecated.
        In general, if you used to call this function, you can now remove that call and still get correct results.
        """
    )
end

# Here charges_h is used as a buffer in the GPU implementation.
# See set_interpolation_points! for details.
# Here `op` is a binary operator `op(new, old)`. For example, to add the new value to a
# previously existent value, pass op = +.
function copy_long_range_output!(
        op::F,
        vs::AbstractVector{<:VectorOfVelocities}, cache::LongRangeCache,
        charges = default_interpolation_output(cache),
    ) where {F}
    (; charges_h,) = cache.common.pointdata
    nout = sum(length, vs)
    nout == length(charges) || throw(DimensionMismatch("wrong length of output vector `vs`"))
    copy_output_values_on_nodes!(op, vs, charges, charges_h)
    vs
end

# These are kept for backwards compatibility for now, but copy_output_values_on_nodes! should be used instead.
function copy_long_range_output!(vs::AbstractVector{<:VectorOfVelocities}, cache::LongRangeCache, args...)
    # By default, only keep the new value, discarding old values in vs.
    op(new, old) = new
    copy_long_range_output!(op, vs, cache, args...)
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

include("backends/exact_sum.jl")
include("backends/nonuniformffts.jl")
