using ..Filaments: FilamentChunkIterator
using ..BiotSavart: BiotSavart, BiotSavartCache, PointData, count_nodes, resize_no_copy!, copy_host_to_device!
using StaticArrays: SVector, SMatrix
using LLVM.Interop: assume
using AcceleratedKernels: AcceleratedKernels as AK
using LinearAlgebra: LinearAlgebra, ⋅, ×

@doc raw"""
    FourierBandForcingBS <: AbstractForcing
    FourierBandForcingBS(; kmin, kmax, α, ε_target, α′ = 0, modify_length = true)

Forcing based on Biot–Savart energetics within a given range of wavenumbers.

This type of forcing does _not_ rely on an imposed normal fluid velocity. Instead, it starts from the
functional derivative of the energy at a given wavevector ``\bm{k}`` with respect to the
vortex positions (according to the Biot–Savart law), and applies a velocity that ensures
energy to increase (if ``α > 0`` or ``ε_{\text{target}} > 0``) at those wavevectors.

One should pass _either_ `α` or `ε_target` but never both. They should be positive for
energy injection (negative values lead to energy dissipation):

- `α` (`\alpha`) is a non-dimensional coefficient which directly sets the amplitude of the forcing velocity;

- `ε_target` (`\varepsilon_target`) has the units of an energy injection rate. In this case,
  the amplitude ``α`` will be adjusted over time in order to keep a roughly constant energy
  injection rate (which in general will _not_ be equal to `ε_target`, see remarks below).

The `modify_length` parameter can be set to `false` to avoid this forcing from increasing
(or decreasing) the vortex length. This changes the forcing velocity by removing the component that is parallel
to the local vortex curvature, such that it doesn't modify the vortex length locally.
This is experimental, and might help (or not) avoiding non-local energy transfers from large
to small scales.

# Extended help

## Forcing definition

The applied forcing velocity has the basic form:

```math
\bm{v}_{\text{f}} = α \bm{s}' × \tilde{\bm{v}}_{\text{s}},
```

where $\tilde{\bm{v}}_{\text{s}}$ is the bandpass filtered Biot–Savart velocity between the
active Fourier wavenumbers `kmin` and `kmax`. This choice ensures that energy is
injected within this range of wavenumbers.

If ``α′ ≠ 0``, an extra term is included which does not modify the energy content at the
chosen wavenumbers:

```math
\bm{v}_{\text{f}} = α \bm{s}' × \tilde{\bm{v}}_{\text{s}} - α′ \bm{s}' × \left[\bm{s}' × \tilde{\bm{v}}_{\text{s}}\right].
```

## Explanations

This forcing attempts to increase the kinetic energy at selected wavenumbers.
Its definition starts from the expression for the kinetic energy at wavenumber ``\bm{k}``:

```math
E(\bm{k}) = \frac{1}{2} |\bm{v}(\bm{k})|^2 = \frac{1}{2k^2} |\bm{ω}(\bm{k})|^2
```

The idea is to translate the vortex positions by ``\bm{s}(ξ) → \bm{s}(ξ) + δ\bm{s}(ξ)`` so
that the energy ``E(\bm{k})`` increases. To determine such a displacement, one can look at
the functional derivative of ``E(\bm{k})`` with respect to the positions ``\bm{s}``:

```math
\begin{align*}
\frac{δE(\bm{k})}{δ\bm{s}}
&= \frac{1}{2k^2} \bm{ω}(\bm{k}) ⋅ \frac{δ\bm{ω}^*(\bm{k})}{δ\bm{s}} + \text{c.c.},
\\
&= \frac{1}{2} \bm{ψ}(\bm{k}) ⋅ \frac{δ\bm{ω}^*(\bm{k})}{δ\bm{s}} + \text{c.c.},
\end{align*}
```

where ``()^*`` denotes a complex conjugate and "c.c." means the complex conjugate of the first term.
Here we have used ``\bm{ω}(\bm{k}) = k^2 \bm{ψ}(\bm{k})``, which corresponds in physical
space to ``\bm{ω} = -∇² \bm{ψ}``.

The functional derivative of vorticity in Fourier space is:

```math
\begin{align*}
\frac{δ\bm{ω}^*(\bm{k})}{δ\bm{s}}
&= \frac{Γ}{V} \frac{δ}{δ\bm{s}} ∮ e^{+i \bm{k} ⋅ \bm{s}(ξ)} \bm{s}(ξ)' \, \text{d}ξ
\\
&= \frac{Γ}{V} \frac{δ}{δ\bm{s}} ∮ \bm{ℒ}[\bm{s}(ξ), \bm{s}(ξ)'] \, \text{d}ξ
\\
&= \frac{Γ}{V} \left[ \frac{∂\bm{ℒ}}{∂\bm{s}} - \frac{\text{d}}{\text{d}ξ} \frac{∂\bm{ℒ}}{∂\bm{s}'} \right]
\\
&= \frac{i Γ}{V} e^{+i \bm{k} ⋅ \bm{s}} \left[ \bm{k} ⊗ \bm{s}' - (\bm{k} ⋅ \bm{s}') I \right] = \frac{i Γ}{V} \bm{B}(\bm{k}, \bm{s}, \bm{s}')
\end{align*}
```

where ``\bm{B}`` is a ``3×3`` complex matrix which has dimensions of an inverse length (``L^{-1}``).

The functional derivative of ``E(\bm{k})`` can now be written as:

```math
\begin{align*}
\frac{δE(\bm{k})}{δ\bm{s}}
&= \frac{Γ}{2V} i \bm{B} \bm{ψ}(\bm{k}) + \text{c.c.}
\\
&= \frac{Γ}{2V} [\bm{s}' × \bm{v}(\bm{k})] e^{i \bm{k} ⋅ \bm{s}} + \text{c.c.}
\\
&= ℜ \left\{ \frac{Γ}{V} [\bm{s}' × \bm{v}(\bm{k})] e^{i \bm{k} ⋅ \bm{s}} \right\}
\end{align*}
```

Finally, the idea is to advect the filaments with a velocity which is parallel to this result for each forced ``\bm{k}``.
One can write this velocity as

```math
\bm{v}_{\text{f}}(\bm{k}, \bm{s}) = α \, ℜ \left\{ [\bm{s}' × \bm{v}(\bm{k})] e^{i \bm{k} ⋅ \bm{s}} \right\} = α \, \bm{v}_0(\bm{k}, \bm{s})
```

where ``α`` is a non-dimensional parameter setting the forcing amplitude.

## Estimating the energy injection rate

One can also try to estimate an energy injection rate at wavevector ``\bm{k}`` associated to this velocity:

```math
\frac{\text{d}E(\bm{k})}{\text{d}t}
= ∮ \frac{δE(\bm{k})}{δ\bm{s}} ⋅ \frac{\text{d}\bm{s}}{\text{d}t} \, \mathrm{d}ξ
= α \frac{Γ}{V} ∮ |\bm{v}_0|^2 \, \mathrm{d}ξ
```

In general, this estimate may be quite inaccurate since the forcing can also affect the
energy at wavevectors other than ``\bm{k}``. But still, this estimate is the one used when
``ε_{\text{target}}`` is given, and may allow to obtain a roughly constant energy injection
rate (even if it's generally different than the "target" one).
"""
struct FourierBandForcingBS{T <: AbstractFloat, N} <: AbstractForcing
    α    :: T
    α′   :: T
    ε_target :: T  # target energy injection rate [L²T⁻³]
    kmin :: T
    kmax :: T
    modify_length :: Bool
    qs   :: Vector{NTuple{N, Int}}  # list of normalised wavevectors (can be empty)
end

function FourierBandForcingBS(; α::Real = 0, ε_target = 0, α′ = 0, modify_length = true, kmin = nothing, kmax = nothing, qs = nothing)
    (α == 0) + (ε_target == 0) == 1 || throw(ArgumentError("one should pass either α or ε_target, but not both"))
    _FourierBandForcingBS(α, α′, ε_target, kmin, kmax, modify_length, qs)
end

function _FourierBandForcingBS(α, α′, ε_target, kmin::Real, kmax::Real, modify_length, qs::Nothing)
    scalars = promote(α, α′, ε_target, kmin, kmax)
    qs = Tuple{}[]  # empty vector of empty tuples
    FourierBandForcingBS(scalars..., modify_length, qs)
end

# This variant allows setting specific q⃗ normalised wavevectors (currently not documented!).
function _FourierBandForcingBS(α, α′, ε_target, kmin::Nothing, kmax::Nothing, modify_length, qs::AbstractVector{NTuple{N, Int}}) where {N}
    kmin = 0
    kmax = 0
    scalars = promote(α, α′, ε_target, kmin, kmax)
    FourierBandForcingBS(scalars..., modify_length, qs)
end

function Base.show(io::IO, f::FourierBandForcingBS{T}) where {T}
    (; α, α′, ε_target, kmin, kmax, qs,) = f
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "FourierBandForcingBS{$T} with:")
    if α != 0
        print(io, "\n$(prefix)├─ Magnitude: α = ", α)
    elseif ε_target != 0
        print(io, "\n$(prefix)├─ Target energy injection rate: ε_target = ", ε_target)
    end
    if α′ != 0
        print(io, "\n$(prefix)├─ Drift coefficient: α′ = ", α′)
    end
    print(io, "\n$(prefix)├─ Forcing modifies length: ", f.modify_length)
    if isempty(qs)
        print(io, "\n$(prefix)└─ Fourier band: |k⃗| ∈ [$kmin, $kmax]")
    else
        print(io, "\n$(prefix)└─ Normalised Fourier wave vectors: |q⃗| = ", qs)
    end
    nothing
end

# Here vs_grid is the superfluid velocity in Fourier space, optionally on a device (GPU).
function init_cache(f::FourierBandForcingBS, cache_bs::BiotSavartCache)
    (; kmin, kmax, qs,) = f
    (; params,) = cache_bs
    vs_grid = BiotSavart.get_longrange_field_fourier(cache_bs).field :: NTuple
    backend = KA.get_backend(vs_grid[1])  # CPU, CUDABackend, ROCBackend, ...
    if isempty(qs)
        # Activate all wavevectors within a Fourier band [kmin, kmax]
        v_h = FourierBandVectorField(undef, params.Ls; kmin, kmax)  # on the host (CPU)
    else
        # Activate selected wavenumbers
        v_h = FourierBandVectorField(undef, params.Ls)  # on the host (CPU)
        for q⃗ in qs
            SyntheticFields.add_normalised_wavevector!(v_h, q⃗)
        end
    end
    v_d = adapt(backend, v_h)  # on the "device" (GPU or CPU)
    prefactor = params.Γ / prod(params.Ls)
    (; v_d, qs = v_d.qs, prefactor,)
end

# After calling this function, cache.v_d contains the Biot-Savart velocity in Fourier space v̂(k⃗).
function _update_cache!(cache, f::FourierBandForcingBS, cache_bs::BiotSavartCache)
    (; v_d,) = cache

    vs_grid, ks_grid, σ_gaussian = let data = BiotSavart.get_longrange_field_fourier(cache_bs)
        local (; state, field, wavenumbers,) = data
        @assert state.quantity == :velocity
        field, wavenumbers, state.smoothing_scale
    end
    @assert σ_gaussian == 0  # fields are unfiltered

    # Copy velocity coefficients within Fourier band
    SyntheticFields.from_fourier_grid!(v_d, vs_grid, ks_grid)

    nothing
end

# This always runs on the CPU.
# This function is adapted from BiotSavart.add_point_charges.
function _compute_geometry!(forcing::FourierBandForcingBS, pointdata_cpu::PointData, fs::AbstractVector{<:ClosedFilament}; quad)
    (; nodes, derivatives_on_nodes, subsegment_lengths) = pointdata_cpu
    Np = count_nodes(fs)
    subsegment_lengths::NTuple{2}
    integration_weights = subsegment_lengths[1]::AbstractVector  # reuse as a buffer
    geom = (; nodes, derivatives_on_nodes, integration_weights)
    # Usually all vectors should already have length Np, but it doesn't harm to try to
    # resize them if that's not the case.
    resize_no_copy!(nodes, Np)
    resize_no_copy!(integration_weights, Np)
    foreach(vs -> resize_no_copy!(vs, Np), derivatives_on_nodes)
    chunks = FilamentChunkIterator(fs)
    @sync for chunk in chunks
        Threads.@spawn for (i, inds, num_nodes_visited) in chunk
            _compute_geometry_filament!(forcing, geom, fs[i], inds, num_nodes_visited; quad)
        end
    end
    geom
end

function _compute_geometry_filament!(::FourierBandForcingBS, geom, f::ClosedFilament, inds::AbstractUnitRange, n::Int; quad)
    (; nodes, derivatives_on_nodes, integration_weights) = geom
    # Obtain length of previous segment assuming closed filament.
    len_prev = let j = first(inds) - 1  # note: we can safely index a filament at i = 0 (same as index i = N)
        sqrt(sum(abs2, f[j + 1] - f[j]))   # length of segment to the right (rough estimate)
    end
    for j in inds
        n += 1
        nodes[n] = f[j]
        derivatives_on_nodes[1][n] = f[j, Derivative(1)]
        derivatives_on_nodes[2][n] = f[j, Derivative(2)]
        len = sqrt(sum(abs2, f[j + 1] - f[j]))    # length of segment to the right (rough estimate)
        # seg = Filaments.Segment(f, j)
        # len = Filaments.segment_length(seg; quad)  # more accurate estimate (but more expensive)
        integration_weights[n] = len_prev + len   # length of the two local segments
        len_prev = len
    end
    n
end

function _to_gpu!(::FourierBandForcingBS, pointdata_gpu::PointData, geom_cpu::NamedTuple)
    (; nodes, derivatives_on_nodes, subsegment_lengths) = pointdata_gpu
    subsegment_lengths::NTuple{2}
    integration_weights = subsegment_lengths[1]::AbstractVector  # reuse as a buffer
    geom_gpu = (; nodes, derivatives_on_nodes, integration_weights)
    foreach(geom_cpu, geom_gpu) do src, dst
        copy_host_to_device!(dst, src)
    end
    geom_gpu
end

function _evaluate_from_geometry!(forcing::FourierBandForcingBS, vf_lin::AbstractVector, geom::NamedTuple, cache)
    (; v_d, prefactor) = cache
    (; nodes, derivatives_on_nodes, integration_weights) = geom
    (; modify_length, ε_target, α, α′) = forcing
    (; qs, cs, Δks,) = v_d

    @assert eachindex(qs) == eachindex(cs)

    resize_no_copy!(vf_lin, length(nodes))

    AK.foreachindex(nodes; block_size = 256) do i
        @inline
        s⃗ = @inbounds nodes[i]
        s⃗_t = @inbounds derivatives_on_nodes[1][i]
        s_t² = sum(abs2, s⃗_t)
        assume(s_t² > 0)
        s_t = sqrt(s_t²)
        s_t_inv = 1 / s_t
        assume(s_t_inv > 0)
        s⃗′ = s⃗_t * s_t_inv  # unit tangent
        V = eltype(vf_lin)  # usually SVector{3, T} == Vec3{T}
        vf = zero(V)  # forcing velocity (excluding α prefactor)
        for n in eachindex(qs)  # iterate over active wavevectors k⃗
            k⃗ = @inbounds SVector(qs[n]) .* Δks
            v̂ = @inbounds cs[n]
            k_dot_s = k⃗ ⋅ s⃗
            local s, c = sincos(k_dot_s)  # note: CUDA defines sincos as well, so it should be fast
            vf_k = real(v̂ * Complex(c, s))
            vf = vf + vf_k
        end
        vf = s⃗′ × vf
        vs_filtered = vf  # velocity filtered at target wavevectors
        if !modify_length
            # If the forcing is not to modify vortex length, we make sure that vf is
            # orthogonal to the local curvature.
            s⃗_tt = @inbounds derivatives_on_nodes[2][i]
            s⃗″ = s⃗_t × (s⃗_tt × s⃗_t)  # curvature vector (up to a scalar constant)
            # @assert s⃗″ ≈ f[j, CurvatureVector()] * (s_t² * s_t²)
            s″_norm2 = sum(abs2, s⃗″)
            assume(s″_norm2 > 0)
            vf = vf - ((vf ⋅ s⃗″) / s″_norm2) * s⃗″
            # @assert norm(vf ⋅ s⃗″) < 1e-12
        end
        if α != 0  # α was prescribed (and possibly α′ as well)
            vf = α * vf - α′ * (s⃗′ × vf)  # this will be the actual forcing velocity
        else
            # Estimate local contribution to energy injection rate.
            # In principle, we should divide dξ by two to get an estimate of the local segment
            # length (taking half the length of each local segment). However, this is
            # compensated by a factor 2 accounting for Hermitian symmetry: if we inject energy into the velocity at
            # a wavevector +k⃗, then we're also injecting the same energy at v(-k⃗).
            dξ = @inbounds integration_weights[i]
            ε_local = (vs_filtered ⋅ vf) * dξ  # estimated energy injection rate around local node (without Γ/V prefactor)
            # We no longer need integration_weights, so we reuse it to store energy injection rates.
            @inbounds integration_weights[i] = ε_local
        end
        @inbounds vf_lin[i] = vf
    end

    if α != 0
        return vf_lin  # nothing else to do, we already have our forcing 
    end

    # If we're here, it's that we prescribed a target ε_inj and we need to adjust the
    # forcing accordingly.
    @assert ε_target != 0

    # Estimate ε_total.
    # Note: integration_weights now contains local contributions to energy injection rate.
    # TODO: pass temporary GPU array to avoid allocations? (see docs for AK.reduce)
    T = eltype(integration_weights)
    ε_total = prefactor * AK.reduce(+, integration_weights; block_size = 256, init = zero(T))

    # If targeting a value of ε, adjust forcing amplitude.
    if ε_total == 0  # just in case; this should never happen in practice
        AK.foreachindex(vf_lin) do i
            @inbounds vf_lin[i] = zero(eltype(vf_lin))
        end
    else
        # In this case, vf_lin only contains the forcing velocities, which need to be
        # rescaled to get the wanted energy injection rate (estimated). We also may need to
        # include the α′ factor if it's nonzero.
        @assert α == 0
        α_actual = ε_target / ε_total
        if α′ == 0
            AK.foreachindex(vf_lin) do i
                @inline
                @inbounds vf_lin[i] = α_actual * vf_lin[i]
            end
        else
            AK.foreachindex(vf_lin) do i
                @inline
                s⃗_t = @inbounds derivatives_on_nodes[1][i]
                s_t² = sum(abs2, s⃗_t)
                assume(s_t² > 0)
                s_t = sqrt(s_t²)
                s_t_inv = 1 / s_t
                assume(s_t_inv > 0)
                s⃗′ = s⃗_t * s_t_inv  # unit tangent
                @inbounds vf_lin[i] = α_actual * vf_lin[i] - α′ * (s⃗′ × vf_lin[i])
            end
        end
    end

    vf_lin
end

function evaluate!(
        forcing::FourierBandForcingBS, cache, vf_all::AbstractVector,
        fs::AbstractVector{<:AbstractFilament}, cache_bs::BiotSavartCache;
        to = TimerOutput(),
    )
    @timeit to "(1) Copy active Fourier coefficients (GPU)" _update_cache!(cache, forcing, cache_bs)

    # To evaluate this forcing, we need:
    #
    # - filament node locations
    # - local tangent vectors
    # - curvature vectors (only when we're targeting an energy injection rate ε_target)
    # - local segment lengths (only when we're targeting an energy injection rate ε_target)
    #
    # We compute them first on the CPU based on filament data, and then we transfer the
    # results to the device where long-range interactions were computed (possibly a GPU),
    # before evaluating the forcing on the GPU. To reduce memory allocations, we reuse
    # temporary arrays (pointdata) defined in the BiotSavart module.
    #
    # TODO: this data may already be in cache_bs.pointdata or cache_bs.longrange.pointdata,
    # so we might be able to skip geometry computation and CPU -> GPU transfers.
    (; quad) = cache_bs.params
    @timeit to "(2) Compute geometry (CPU)" geom_cpu = _compute_geometry!(forcing, cache_bs.pointdata, fs; quad)::NamedTuple
    @timeit to "(3) Copy geometry CPU -> GPU" geom_gpu = _to_gpu!(forcing, cache_bs.longrange.pointdata, geom_cpu)::NamedTuple

    # Forcing velocities in "linear" format.
    vf_lin = BiotSavart.default_interpolation_output(cache_bs.longrange)  # reused as temporary array to hold forcing velocity
    @timeit to "(4) Compute forcing (GPU)" _evaluate_from_geometry!(forcing, vf_lin, geom_gpu, cache)

    @timeit to "(5) Copy forcing GPU -> CPU" begin
        vs_buf = cache_bs.pointdata.nodes  # CPU array used as intermediate buffer
        BiotSavart.copy_output_values_on_nodes!(vf_all, vf_lin, vs_buf)
    end

    vf_all
end
