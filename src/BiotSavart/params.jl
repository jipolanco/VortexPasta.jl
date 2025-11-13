using ..Constants: RealConst

const MaybeConst{T} = Union{RealConst, T}

# Common parameters to short- and long-range computations.
# Note: it would be nice to further constrain the Alpha parameter (commented code below),
# but then precompilation fails (enters an infinite loop?). Tested on Julia 1.10-rc1.
struct ParamsCommon{
        T <: AbstractFloat,
        Alpha <: Real,  #  <: MaybeConst{T} (fails!!)
        Sigma <: Real,  #  <: MaybeConst{T} (fails!!)
        Periods <: NTuple{3, MaybeConst{T}},
        Quad <: StaticSizeQuadrature,
        QuadMaybeAdaptive <: Union{StaticSizeQuadrature, PreallocatedQuadrature{T}},
    }
    Γ  :: T        # vortex circulation
    a  :: T        # vortex core size
    Δ  :: T        # LIA coefficient given by core vorticity profile
    α  :: Alpha    # Ewald splitting parameter (inverse length scale)
    σ  :: Sigma    # Ewald splitting length scale = 1 / α√2 = std of Gaussian filter
    Ls :: Periods  # size of unit cell (= period in each direction)
    quad :: Quad   # quadrature rule used for short- and long-range computations
    quad_near_singularity :: QuadMaybeAdaptive  # quadrature rule to be used near singularities (adaptive by default)
    avoid_explicit_erf :: Bool
    function ParamsCommon{T}(Γ, a, Δ, α_in, Ls_in, quad, quad_near_singularity, avoid_explicit_erf) where {T}
        α = maybe_convert(T, α_in)  # don't convert constants (e.g. if α = Zero())
        Ls = map(L -> maybe_convert(T, L), Ls_in)
        σ = 1 / (α * sqrt(T(2)))
        quad_ns = convert(T, quad_near_singularity)  # does nothing if the quadrature has the right type
        new{T, typeof(α), typeof(σ), typeof(Ls), typeof(quad), typeof(quad_ns)}(
            Γ, a, Δ, α, σ, Ls, quad, quad_ns, avoid_explicit_erf,
        )
    end
end

maybe_convert(::Type{T}, x::Real) where {T <: AbstractFloat} = convert(T, x)
maybe_convert(::Type{T}, x::RealConst) where {T <: AbstractFloat} = x  # don't convert constants (Zero, Infinity)

function Base.show(io::IO, p::ParamsCommon)
    (; Γ, a, Δ, Ls, α,) = p
    σ = 1 / (α * sqrt(2))
    print(io, "\n - Physical parameters:")
    print(io, "\n   * Vortex circulation:         Γ  = ", Γ)
    print(io, "\n   * Vortex core radius:         a  = ", a)
    print(io, "\n   * Vortex core parameter:      Δ  = ", Δ)
    print(io, "\n   * Domain period:              Ls = ", Ls)
    print(io, "\n - Numerical parameters:")
    print(io, "\n   * Ewald splitting parameter:  α = ", α, " (σ = 1/α√2 = ", σ, ")")
    print(io, "\n   * Quadrature rule:            ", p.quad)
    print(io, "\n   * Quadrature rule (alt.):     ", p.quad_near_singularity, " (used near singularities)")
    print(io, "\n   * Avoid explicit erf:         ", p.avoid_explicit_erf)
    nothing
end

function to_hdf5(g, p::ParamsCommon{T}) where {T}
    g["Gamma"] = p.Γ
    g["a"] = p.a
    g["Delta"] = p.Δ
    g["Ls"] = collect(T, p.Ls)    # this allows converting Infinity() to a float (e.g. Inf)
    g["alpha"] = convert(T, p.α)  # convert Zero() -> 0.0
    g["quad"] = string(p.quad)
    g["quad_near_singularity"] = string(p.quad_near_singularity)
    g["avoid_explicit_erf"] = string(p.avoid_explicit_erf)
    nothing
end

Base.eltype(::Type{<:ParamsCommon{T}}) where {T} = T
Base.eltype(p::ParamsCommon) = eltype(typeof(p))

function Base.convert(::Type{T}, p::ParamsCommon) where {T <: AbstractFloat}
    (; Γ, a, Δ, α, Ls, quad, quad_near_singularity, avoid_explicit_erf) = p
    ParamsCommon{T}(Γ, a, Δ, α, Ls, quad, quad_near_singularity, avoid_explicit_erf)  # converts all floats to type T
end

## ================================================================================ ##

struct ParamsShortRange{
        T <: AbstractFloat,
        Backend <: ShortRangeBackend,
        Quadrature <: StaticSizeQuadrature,
        Common <: ParamsCommon{T},
        CutoffDist <: MaybeConst{T},
        LIASegmentFraction <: Union{Nothing, T},
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    rcut    :: CutoffDist  # cutoff distance
    rcut_sq :: CutoffDist  # squared cutoff distance
    lia_segment_fraction :: LIASegmentFraction
    use_simd :: Bool

    function ParamsShortRange(
            backend::ShortRangeBackend, quad::StaticSizeQuadrature,
            common::ParamsCommon{T}, rcut_::Real,
            lia_segment_fraction_in,
            use_simd::Bool,
        ) where {T}
        (; Ls,) = common
        rcut = maybe_convert(T, rcut_)
        lia_segment_fraction = lia_segment_fraction_in === nothing ? nothing : convert(T, lia_segment_fraction_in)
        rcut_max = max_cutoff_distance(backend, Ls)
        rcut > rcut_max && error(
            lazy"""cutoff distance `rcut = $rcut` is larger than that allowed by the $(typeof(backend)) backend.
            See docs for $(typeof(backend)) for more details.
            """,
        )
        2 * rcut ≤ min(Ls...) ||
            error(lazy"cutoff distance `rcut = $rcut` is too large. It must be less than half the cell unit size `L` in each direction: Ls = $Ls.")
        rcut_sq = rcut * rcut
        new{
            T, typeof(backend), typeof(quad), typeof(common), typeof(rcut),
            typeof(lia_segment_fraction)
        }(
            backend, quad, common, rcut, rcut_sq, lia_segment_fraction, use_simd,
        )
    end
end

function Base.show(io::IO, p::ParamsShortRange)
    (; common, rcut, lia_segment_fraction,) = p
    (; Ls, α,) = common
    β_shortrange = α === Zero() ? (rcut === Infinity() ? Infinity() : Zero()) : α * rcut
    rcut_L = rcut === Infinity() ? rcut : rcut / minimum(Ls)  # avoid Infinity() / Infinity()
    print(io, "\n   * Short-range backend:        ", p.backend)
    print(io, "\n   * Short-range cut-off:        r_cut = ", rcut, " (r_cut/L = ", rcut_L, ")")
    print(io, "\n   * Short-range cut-off coeff.: β_shortrange = ", β_shortrange)
    print(io, "\n   * Local segment fraction:     ", something(lia_segment_fraction, 1))
    print(io, "\n   * Short-range uses SIMD:      ", p.use_simd)
    nothing
end

function to_hdf5(g, p::ParamsShortRange{T}) where {T}
    (; common, rcut, lia_segment_fraction,) = p
    (; α,) = common
    β_shortrange = α === Zero() ? (rcut === Infinity() ? Infinity() : Zero()) : α * rcut
    g["backend_shortrange"] = string(p.backend)
    g["beta_shortrange"] = β_shortrange
    g["r_cut"] = convert(T, rcut)
    g["lia_segment_fraction"] = lia_segment_fraction === nothing ? one(T) : lia_segment_fraction  # setting this to nothing is equivalent to setting it to 1
    g["use_simd"] = p.use_simd
    nothing
end

# Create a new ParamsShortRange based on the type of ParamsCommon.
function change_float_type(p::ParamsShortRange, common::ParamsCommon{T}) where {T}
    ParamsShortRange(p.backend, p.quad, common, p.rcut, p.lia_segment_fraction, p.use_simd)
end

backend(p::ParamsShortRange) = p.backend
quadrature_rule(p::ParamsShortRange) = p.quad

## ================================================================================ ##

struct ParamsLongRange{
        T <: AbstractFloat,
        Backend <: LongRangeBackend,
        Quadrature <: StaticSizeQuadrature,
        Common <: ParamsCommon{T},
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    Ns      :: Dims{3}     # grid dimensions for FFTs
    truncate_spherical :: Bool  # if true, perform spherical truncation in Fourier space
end

has_real_to_complex(p::ParamsLongRange) = has_real_to_complex(p.backend)

function maximum_wavenumber(p::ParamsLongRange{T}) where {T}
    minimum(zip(p.Ns, p.common.Ls)) do (N, L)
        local m = (N - 1) ÷ 2
        convert(T, 2π * m / L)  # convert in case m/L = Zero()
    end
end

function Base.show(io::IO, p::ParamsLongRange{T}) where {T}
    (; common, Ns,) = p
    (; α,) = common
    kmax = maximum_wavenumber(p)
    β_longrange = kmax === Zero() ? Zero() : kmax / (2 * α)
    print(io, "\n   * Long-range backend:         ", p.backend)
    print(io, "\n   * Long-range resolution:      Ns = ", Ns, " (kmax = ", kmax, ")")
    print(io, "\n   * Long-range cut-off coeff.:  β_longrange = ", β_longrange)
    print(io, "\n   * Long-range spherical truncation: ", p.truncate_spherical)
    nothing
end

function to_hdf5(g, p::ParamsLongRange{T}) where {T}
    (; common, Ns,) = p
    (; α,) = common
    kmax = maximum_wavenumber(p)
    β_longrange = kmax === Zero() ? Zero() : kmax / (2 * α)
    g["backend_longrange"] = string(p.backend)
    g["beta_longrange"] = β_longrange
    g["Ns"] = collect(Ns)
    g["truncate_spherical"] = p.truncate_spherical
    nothing
end

function change_float_type(p::ParamsLongRange, common::ParamsCommon{T}) where {T}
    ParamsLongRange(p.backend, p.quad, common, p.Ns, p.truncate_spherical)
end

backend(p::ParamsLongRange) = p.backend
quadrature_rule(p::ParamsLongRange) = p.quad

## ================================================================================ ##

"""
    ParamsBiotSavart{T <: AbstractFloat}

Contains parameters for calculation of Biot–Savart integrals using fast Ewald splitting.

The type parameter `T` corresponds to the precision used in computations
(typically `Float64` or `Float32`).

# Construction

    ParamsBiotSavart([T = Float64]; Γ, a, α, Ls, Ns, rcut, optional_kws...)

where the optional parameter `T` sets the numerical precision.

Mandatory and optional keyword arguments are detailed in the extended help below.

See also [`BiotSavart.autotune`](@ref) for an alternative way of setting Biot–Savart parameters.

# Extended help

## Mandatory keyword arguments

- `Γ::Real`: vortex circulation (assumed constant);

- `a::Real`: vortex core size (assumed constant);

- `α::Real`: Ewald splitting parameter (inverse length scale). One can set
  `α = Zero()` to efficiently disable long-range computations.

- `Ls::Union{Real, NTuple{3, Real}}`: domain period in each Cartesian direction.
  If a single value is passed (e.g. `Ls = 2π`), it is assumed that periods are
  the same in all directions.

  One can set `Ls = Infinity()` to disable periodicity. This should be done in combination with `α = Zero()`.

- `Ns::Dims{3}`: dimensions of physical grid used for long-range interactions.
  This parameter is not required if `α = Zero()`.

- `rcut`: cutoff distance for computation of short-range interactions.
  For performance and practical reasons, the cutoff distance must be less than half the cell
  unit size in each direction, i.e. `rcut < minimum(Ls) / 2`.
  This parameter is not required if `α = Zero()`.

## Optional keyword arguments and their defaults

### General

- `quadrature::StaticSizeQuadrature = GaussLegendre(3)`: quadrature rule for short- and
  long-range interactions. For example, if `quadrature = GaussLegendre(4)`, then 4 evaluations
  of the Biot–Savart integrand will be done for each filament segment.

### Short-range interactions

- `backend_short::ShortRangeBackend`: backend used to compute
  short-range interactions. The default is `CellListsBackend(2)`, unless periodicity is
  disabled, in which case `NaiveShortRangeBackend()` is used.
  See [`ShortRangeBackend`](@ref) for a list of possible backends.

### Long-range interactions

- `backend_long::LongRangeBackend = NonuniformFFTsBackend()`: backend used to compute
  long-range interactions. See [`LongRangeBackend`](@ref) for a list of possible backends;

- `longrange_truncate_spherical = false`: if `true`, perform a spherical truncation in
  Fourier space, discarding all wavenumbers such that ``|\\bm{k}| > k_{\\text{max}}``.
  This is not recommended as it leads to precision loss, and should be used for testing only
  (in particular, for verifying error estimates which assume this kind of truncation).

### Local self-induced velocity

- `Δ = 0.25`: coefficient appearing in the local self-induced velocity (LIA
  term), which depends on the vorticity profile at the vortex core.

  Some common values of `Δ` are:

  * `Δ = 0.25` for a constant vorticity profile (default);

  * `Δ = 0.5` for a hollow vortex;

  * `Δ ≈ 0.905 ≈ 0.558 + ln(2) / 2` for a Gaussian vorticity profile (taking
    `a` as the Gaussian standard deviation `σ`);

  * `Δ ≈ 0.615` for a Gross–Pitaevskii vortex with healing length `a`.

  See [Saffman1993](@citet), sections 10.2--10.3 for the first three.

- `lia_segment_fraction = nothing`: can be used to indicate that the LIA term should be
  evaluated over a *fraction* of the two segments surrounding a node. In this case, it
  should be a real value in ``(0, 1]``. The default (`nothing`) is equivalent to 1, and
  means that the LIA term is evaluated over the full segments. If smaller than 1, the
  velocity induced by the excluded part of the segments will be evaluated using the regular
  Biot–Savart law (using quadratures within each subsegment). This may improve accuracy,
  especially when the discretisation distance is relatively large. See also
  `quadrature_near_singularity` for the choice of quadrature rule when this is enabled.

- `quadrature_near_singularity = quadrature`: quadrature rule to be used
  when integrating near a singularity when `lia_segment_fraction` is enabled.
  By default this is equal to the quadrature rule used for non-local interactions (the
  `quadrature` parameter), which has static size and is not adaptive. For even higher
  accuracy, one may choose an adaptive quadrature rule such as [`AdaptiveTanhSinh`](@ref),
  e.g. `quadrature_near_singularity = AdaptiveTanhSinh(T; nlevels = 5)`, but note that this
  is costly and might not lead to considerable accuracy gains.

### Other performance parameters

- `use_simd::Bool = true`: whether to use explicit SIMD during the computation of short-range
  interactions. This applies to the CPU implementation of short-range interactions in
  periodic domains only, where the SIMD implementation can accelerate the computation of
  `erfc(αr)` in particular.
  Usually there is no reason to disable this, other than to verify the accuracy or
  performance of the SIMD implementation;

- `avoid_explicit_erf::Bool = true`: whether to avoid explicit computation of `erf(αr)`,
  which appears when one wants to subtract the local contribution of long-range interactions
  (the "self"-interaction), and which may be a bit expensive. If this is `true`, this is
  avoided by _including_ the local segments in short-range interactions, and then
  subtracting the _full_ Biot-Savart integrals (instead of long-range ones) evaluated over
  the local segments. Note that the local integrals (short-range and full) are mathematically singular,
  and numerically might lead to some loss of accuracy due to [catastrophic cancellation](https://en.wikipedia.org/wiki/Catastrophic_cancellation),
  especially when working in single precision (`T = Float32`).
"""
struct ParamsBiotSavart{
        T,
        Common <: ParamsCommon{T},
        ShortRange <: ParamsShortRange,
        LongRange <: ParamsLongRange,
    }
    common     :: Common
    shortrange :: ShortRange
    longrange  :: LongRange
end

function ParamsBiotSavart(
        ::Type{T}, Γ::Real, α::Real, Ls::NTuple{3, Real};
        a::Real,
        quadrature::StaticSizeQuadrature = GaussLegendre(3),
        quadrature_near_singularity::AbstractQuadrature = quadrature,
        backend_short::ShortRangeBackend = default_short_range_backend(Ls),
        backend_long::LongRangeBackend = default_long_range_backend(Ls),
        longrange_truncate_spherical::Bool = false,
        Δ::Real = 0.25,
        lia_segment_fraction::Union{Nothing, Real} = nothing,
        use_simd::Bool = default_use_simd(Ls),
        avoid_explicit_erf::Bool = true,
        kws...,
    ) where {T}
    # TODO better split into physical (Γ, a, Δ, Ls) and numerical (α, rcut, Ns, ...) parameters?
    # - define ParamsPhysical instead of ParamsCommon
    # - include α in both ParamsShortRange and ParamsLongRange?
    (; Ns, rcut,) = _extra_params(α; kws...)
    if lia_segment_fraction !== nothing && quadrature_near_singularity === NoQuadrature()
        @warn(
            "quadrature_near_singularity has been set to NoQuadrature(); this will lead to wrong results!",
            quadrature, quadrature_near_singularity, lia_segment_fraction  # print these variables for extra information
        )
    end
    common = ParamsCommon{T}(Γ, a, Δ, α, Ls, quadrature, quadrature_near_singularity, avoid_explicit_erf)
    sr = ParamsShortRange(backend_short, quadrature, common, rcut, lia_segment_fraction, use_simd)
    lr = ParamsLongRange(backend_long, quadrature, common, Ns, longrange_truncate_spherical)
    ParamsBiotSavart(common, sr, lr)
end

function Base.convert(::Type{T}, p::ParamsBiotSavart) where {T <: AbstractFloat}
    common = convert(T, p.common)
    shortrange = change_float_type(p.shortrange, common)
    longrange = change_float_type(p.longrange, common)
    ParamsBiotSavart(common, shortrange, longrange)
end

# Returns `true` if `Ls` contains `Infinity` (one or more times), `false` otherwise.
is_open_domain(Ls::Tuple) = is_open_domain(Ls...)
is_open_domain(::Infinity, etc...) = true
is_open_domain(::Real, etc...) = is_open_domain(etc...)
is_open_domain() = false

function default_short_range_backend(Ls::Tuple)
    if is_open_domain(Ls)
        NaiveShortRangeBackend()
    else
        CellListsBackend(2)  # use 2 subdivisions by default (generally faster)
    end
end

function default_long_range_backend(Ls::Tuple)
    if is_open_domain(Ls)
        NullLongRangeBackend()
    else
        NonuniformFFTsBackend(σ = 1.5, m = HalfSupport(4))
    end
end

# By default, use explicit SIMD in periodic domains (this option is simply ignored in
# non-periodic domains, as explicit SIMD is not implemented).
default_use_simd(Ls::Tuple) = is_open_domain(Ls) ? false : true

# Returns the float type used (e.g. Float64)
Base.eltype(::Type{<:ParamsBiotSavart{T}}) where {T} = T
Base.eltype(p::ParamsBiotSavart) = eltype(typeof(p))

"""
    BiotSavart.circulation(p::ParamsBiotSavart) -> Γ

Return the circulation `Γ` associated to each vortex.
"""
circulation(p::ParamsBiotSavart) = p.common.Γ

"""
    BiotSavart.periods(p::ParamsBiotSavart) -> (Lx, Ly, Lz)

Return the domain periods in each direction.
"""
periods(p::ParamsBiotSavart) = p.common.Ls

"""
    BiotSavart.domain_is_periodic(p::ParamsBiotSavart) -> Bool
    BiotSavart.domain_is_periodic(Ls::NTuple) -> Bool

Check whether the domain is periodic.

Returns `true` if the domain is periodic in *all* directions, `false` otherwise.
"""
domain_is_periodic(p::ParamsBiotSavart) = domain_is_periodic(periods(p))
domain_is_periodic(Ls::NTuple) = !is_open_domain(Ls)

_extra_params(α::Zero; Ns = nothing, rcut = ∞) = (; Ns = (0, 0, 0), rcut,)  # Ns is always (0, 0, 0), no matter the input
_extra_params(α::Real; Ns, rcut,) = (; Ns, rcut,)  # Ns and rcut are required in this case

ParamsBiotSavart(::Type{T}; Γ::Real, α::Real, Ls, kws...) where {T} =
    ParamsBiotSavart(T, Γ, α, _convert_periods(Ls); kws...)

_convert_periods(Ls::NTuple{3, Real}) = Ls
_convert_periods(L::Real) = (L, L, L)

ParamsBiotSavart(; kws...) = ParamsBiotSavart(Float64; kws...)

# This is for convenience: doing p.α is equivalent to p.common.α.
@inline function Base.getproperty(p::ParamsBiotSavart, name::Symbol)
    common = getfield(p, :common)
    if hasproperty(common, name)
        getproperty(common, name)
    else
        getfield(p, name)
    end
end

function Base.propertynames(p::ParamsBiotSavart, private::Bool = false)
    (fieldnames(typeof(p))..., propertynames(p.common, private)...)
end

function Base.show(io::IO, p::ParamsBiotSavart{T}) where {T}
    print(io, "ParamsBiotSavart{$T} with:")
    show(io, p.common)
    show(io, p.shortrange)
    show(io, p.longrange)
    nothing
end

# Write parameters to HDF5 output (or anything that behaves like an HDF5 file/group.
function to_hdf5(g, p::ParamsBiotSavart)
    to_hdf5(g, p.common)
    to_hdf5(g, p.shortrange)
    to_hdf5(g, p.longrange)
    nothing
end

function Base.summary(io::IO, p::ParamsBiotSavart{T}) where {T}
    # Print a few physical parameters (and α)
    (; Γ, a, Δ, α,) = p.common
    print(io, "ParamsBiotSavart{$T}(Γ = $Γ, a = $a, Δ = $Δ, α = $α, …)")
end

@doc raw"""
    BiotSavart.kelvin_wave_period(p::ParamsBiotSavart, λ::Real) -> Real

Return the period ``T(λ)`` associated to Kelvin waves of wavelength ``λ``.

This can be convenient for setting the timestep `dt` associated to a filament discretisation
distance `δ`. The timestep should typically be proportional to the period of the Kelvin
waves of wavelength `δ`.

The Kelvin wave period is ``T(λ) = 2π/ω(k)`` where ``k = 2π/λ`` is the wavenumber associated
to ``λ`` and ``ω(k)`` is the Kelvin wave dispersion relation:

```math
ω(k) = \frac{Γ k^2}{4π} \left[
  \ln\left( \frac{2}{k a} \right) - γ + \frac{1}{2} - Δ
\right]
```

where ``γ ≈ 0.5772`` is the Euler--Mascheroni constant.
"""
function kelvin_wave_period end

# Returns the expected period of a small-amplitude Kelvin wave of wavelength λ.
# Can be useful when setting a simulation timestep.
kelvin_wave_period(λ::Real; a, Δ, Γ) = 2 * λ^2 / Γ / (
    log(λ / (π * a)) + 1/2 - (Δ + MathConstants.γ)
)
kelvin_wave_period(p::ParamsBiotSavart, λ::Real) = kelvin_wave_period(λ; a = p.a, Δ = p.Δ, Γ = p.Γ)
