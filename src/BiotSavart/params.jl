using ..BasicTypes: RealConst

const MaybeConst{T} = Union{RealConst, T}

# Common parameters to short- and long-range computations.
# Note: it would be nice to further constrain the Alpha parameter (commented code below),
# but then precompilation fails (enters an infinite loop?). Tested on Julia 1.10-rc1.
struct ParamsCommon{
        T <: AbstractFloat,
        Alpha <: Real,  #  <: MaybeConst{T} (fails!!)
        Sigma <: Real,  #  <: MaybeConst{T} (fails!!)
        Periods <: NTuple{3, MaybeConst{T}},
        Quad <: AbstractQuadrature,
    }
    Γ  :: T        # vortex circulation
    a  :: T        # vortex core size
    Δ  :: T        # LIA coefficient given by core vorticity profile
    α  :: Alpha    # Ewald splitting parameter (inverse length scale)
    σ  :: Sigma    # Ewald splitting length scale = 1 / α√2 = std of Gaussian filter
    Ls :: Periods  # size of unit cell (= period in each direction)
    quad :: Quad   # quadrature rule used for short- and long-range computations
    function ParamsCommon{T}(Γ, a, Δ, α_in, Ls_in, quad) where {T}
        α = maybe_convert(T, α_in)  # don't convert constants (e.g. if α = Zero())
        Ls = map(L -> maybe_convert(T, L), Ls_in)
        σ = 1 / (α * sqrt(T(2)))
        new{T, typeof(α), typeof(σ), typeof(Ls), typeof(quad)}(Γ, a, Δ, α, σ, Ls, quad)
    end
end

maybe_convert(::Type{T}, x::Real) where {T <: AbstractFloat} = convert(T, x)
maybe_convert(::Type{T}, x::RealConst) where {T <: AbstractFloat} = x  # don't convert constants (Zero, Infinity)

function Base.show(io::IO, p::ParamsCommon)
    (; Γ, a, Δ, Ls, α, quad,) = p
    σ = 1 / (α * sqrt(2))
    print(io, "\n - Physical parameters:")
    print(io, "\n   * Vortex circulation:         Γ  = ", Γ)
    print(io, "\n   * Vortex core radius:         a  = ", a)
    print(io, "\n   * Vortex core parameter:      Δ  = ", Δ)
    print(io, "\n   * Domain period:              Ls = ", Ls)
    print(io, "\n - Numerical parameters:")
    print(io, "\n   * Ewald splitting parameter:  α = ", α, " (σ = 1/α√2 = ", σ, ")")
    print(io, "\n   * Quadrature rule:            ", quad)
    nothing
end

Base.eltype(::Type{<:ParamsCommon{T}}) where {T} = T
Base.eltype(p::ParamsCommon) = eltype(typeof(p))

function Base.convert(::Type{T}, p::ParamsCommon) where {T <: AbstractFloat}
    (; Γ, a, Δ, α, Ls, quad,) = p
    ParamsCommon{T}(Γ, a, Δ, α, Ls, quad)  # converts all floats to type T
end

## ================================================================================ ##

struct ParamsShortRange{
        T <: Real,
        Backend <: ShortRangeBackend,
        Quadrature <: AbstractQuadrature,
        Common <: ParamsCommon{T},
        CutoffDist <: MaybeConst{T},
        RegulariseBinormal <: StaticBool,
        LIASegmentFraction <: Union{Nothing, Real}
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    rcut    :: CutoffDist  # cutoff distance
    rcut_sq :: CutoffDist  # squared cutoff distance
    regularise_binormal  :: RegulariseBinormal
    lia_segment_fraction :: LIASegmentFraction

    function ParamsShortRange(
            backend::ShortRangeBackend, quad::AbstractQuadrature,
            common::ParamsCommon{T}, rcut_::Real, regularise::StaticBool,
            lia_segment_fraction,
        ) where {T}
        (; Ls,) = common
        rcut = maybe_convert(T, rcut_)
        2 * rcut ≤ min(Ls...) ||
            error(lazy"cutoff distance `rcut = $rcut` is too large. It must be less than half the cell unit size `L` in each direction: Ls = $Ls.")
        rcut_sq = rcut * rcut
        new{
            T, typeof(backend), typeof(quad), typeof(common), typeof(rcut),
            typeof(regularise), typeof(lia_segment_fraction)
        }(
            backend, quad, common, rcut, rcut_sq, regularise, lia_segment_fraction,
        )
    end
end

function Base.show(io::IO, p::ParamsShortRange)
    (; common, rcut,) = p
    (; Ls, α,) = common
    β_shortrange = α === Zero() ? (rcut === Infinity() ? Infinity() : Zero()) : α * rcut
    rcut_L = rcut === Infinity() ? rcut : rcut / minimum(Ls)  # avoid Infinity() / Infinity()
    print(io, "\n   * Short-range backend:        ", p.backend)
    print(io, "\n   * Short-range cut-off:        r_cut = ", rcut, " (r_cut/L = ", rcut_L, ")")
    print(io, "\n   * Short-range cut-off coeff.: β_shortrange = ", β_shortrange)
    nothing
end

# Create a new ParamsShortRange based on the type of ParamsCommon.
function change_float_type(p::ParamsShortRange, common::ParamsCommon{T}) where {T}
    ParamsShortRange(p.backend, p.quad, common, p.rcut, p.regularise_binormal, p.lia_segment_fraction)
end

struct ParamsLongRange{
        T,
        Backend <: LongRangeBackend,
        Quadrature <: AbstractQuadrature,
        Common <: ParamsCommon{T},
    }
    backend :: Backend
    quad    :: Quadrature  # quadrature rule used for numerical integration
    common  :: Common      # common parameters (Γ, α, Ls)
    Ns      :: Dims{3}     # grid dimensions for FFTs
end

function Base.show(io::IO, p::ParamsLongRange{T}) where {T}
    (; common, Ns,) = p
    (; α, Ls,) = common
    kmax = minimum(zip(Ns, Ls)) do (N, L)
        local m = (N - 1) ÷ 2
        convert(T, 2π * m / L)  # convert in case m/L = Zero()
    end
    β_longrange = kmax === Zero() ? Zero() : kmax / (2 * α)
    print(io, "\n   * Long-range backend:         ", p.backend)
    print(io, "\n   * Long-range resolution:      Ns = ", Ns, " (kmax = ", kmax, ")")
    print(io, "\n   * Long-range cut-off coeff.:  β_longrange = ", β_longrange)
    nothing
end

function change_float_type(p::ParamsLongRange, common::ParamsCommon{T}) where {T}
    ParamsLongRange(p.backend, p.quad, common, p.Ns)
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

    ParamsBiotSavart([T = Float64]; Γ, a, α, Ls, Ns, optional_kws...)

where the optional parameter `T` sets the numerical precision.

Mandatory and optional keyword arguments are detailed in the extended help below.

# Extended help

## Mandatory keyword arguments

- `Γ::Real`: vortex circulation (assumed constant);

- `a::Real`: vortex core size (assumed constant);

- `α::Real`: Ewald splitting parameter (inverse length scale). One can set
  `α = Zero()` to efficiently disable long-range computations.

- `Ls::Union{Real, NTuple{3, Real}}`: size of unit cell (i.e. period in each direction).
  If a single value is passed (e.g. `Ls = 2π`), it is assumed that periods are
  the same in each direction.

  One can set `Ls = Infinity()` to disable periodicity. This should be done in combination with `α = Zero()`.

- `Ns::Dims{3}`: dimensions of physical grid used for long-range interactions. This parameter
  is not required if `α = Zero()`.

## Optional keyword arguments (and their defaults)

### General

- `quadrature::AbstractQuadrature = GaussLegendre(3)`: quadrature rule for short- and
  long-range interactions. For example, if `quadrature = GaussLegendre(4)`, then 4 evaluations
  of the Biot–Savart integrand will be done for each filament segment.

### Short-range interactions

- `backend_short::ShortRangeBackend`: backend used to compute
  short-range interactions. The default is `CellListsBackend(2)`, unless periodicity is
  disabled, in which case `NaiveShortRangeBackend()` is used.
  See [`ShortRangeBackend`](@ref) for a list of possible backends;

- `rcut = 4√2 / α`: cutoff distance for computation of short-range interactions.
  For performance and practical reasons, the cutoff distance must be less than half the cell
  unit size in each direction, i.e. `rcut < minimum(Ls) / 2`.

### Long-range interactions

- `backend_long::LongRangeBackend = NonuniformFFTsBackend()`: backend used to compute
  long-range interactions. See [`LongRangeBackend`](@ref) for a list of possible backends.

### Local self-induced velocity

- `Δ = 0.25`: coefficient appearing in the local self-induced velocity (LIA
  term), which depends on the vorticity profile at the vortex core.

  Some common values of `Δ` are:

  * `Δ = 0.25` for a constant vorticity profile (default);

  * `Δ = 0.5` for a hollow vortex;

  * `Δ ≈ 0.905 ≈ 0.558 + ln(2) / 2` for a Gaussian vorticity profile (taking
    `a` as the Gaussian standard deviation `σ`);

  * `Δ ≈ 0.615` for a Gross–Pitaevskii vortex with healing length `a`.

  See Saffman (1992), sections 10.2--10.3 for the first three.

- `lia_segment_fraction = nothing`: can be used to indicate that the LIA term should be
  evaluated over a *fraction* of the two segments surrounding a node. In this case, it
  should be a real value in ``(0, 1]``. The default (`nothing`) is equivalent to 1, and
  means that the LIA term is evaluated over the full segments. If smaller than 1, the
  velocity induced by the excluded part of the segments will be evaluated using the regular
  Biot–Savart law (using quadratures within each subsegment). This may improve accuracy,
  especially when the discretisation distance is relatively large.

- `regularise_binormal = Val(false)`: if `Val(true)`, regularise the estimation of the local
  binormal vector for computation of the LIA term. The binormal vector is averaged along the
  local filament segments. This may lead to more stable simulations of single vortex rings,
  but it's not recommended (nor very useful) for more complex cases. In particular, it can
  lead to spurious energy fluctuations.

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
        quadrature::AbstractQuadrature = GaussLegendre(3),
        quadrature_short = nothing,  # deprecated
        quadrature_long = nothing,   # deprecated
        backend_short::ShortRangeBackend = default_short_range_backend(Ls),
        backend_long::LongRangeBackend = NonuniformFFTsBackend(),
        regularise_binormal::Val{RegulariseBinormal} = Val(false),
        Δ::Real = 0.25,
        lia_segment_fraction::Union{Nothing, Real} = nothing,
        kws...,
    ) where {T, RegulariseBinormal}
    # TODO better split into physical (Γ, a, Δ, Ls) and numerical (α, rcut, Ns, ...) parameters?
    # - define ParamsPhysical instead of ParamsCommon
    # - include α in both ParamsShortRange and ParamsLongRange?
    (; Ns, rcut,) = _extra_params(α; kws...)
    if !isnothing(quadrature_short) || !isnothing(quadrature_long)
        @warn "`quadrature_short` and `quadrature_long` are deprecated and will be removed. Pass `quadrature` instead."
    end
    quad = _parse_quadrature_args(quadrature, quadrature_short, quadrature_long)
    common = ParamsCommon{T}(Γ, a, Δ, α, Ls, quad)
    sr = ParamsShortRange(
        backend_short, quad, common, rcut, StaticBool(RegulariseBinormal),
        lia_segment_fraction,
    )
    lr = ParamsLongRange(backend_long, quad, common, Ns)
    ParamsBiotSavart(common, sr, lr)
end

function Base.convert(::Type{T}, p::ParamsBiotSavart) where {T <: AbstractFloat}
    common = convert(T, p.common)
    shortrange = change_float_type(p.shortrange, common)
    longrange = change_float_type(p.longrange, common)
    ParamsBiotSavart(common, shortrange, longrange)
end

_parse_quadrature_args(quad::AbstractQuadrature, ::Nothing, ::Nothing) = quad
_parse_quadrature_args(::AbstractQuadrature, short::AbstractQuadrature, ::Nothing) = short
_parse_quadrature_args(::AbstractQuadrature, ::Nothing, long::AbstractQuadrature) = long
function _parse_quadrature_args(::AbstractQuadrature, short::AbstractQuadrature, long::AbstractQuadrature)
    short === long || throw(ArgumentError("`quadrature_short` and `quadrature_long` must be equal. Use `quadrature` instead to set both."))
    short
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

Check whether the domain is periodic.

Returns `true` if the domain is periodic in *all* directions, `false` otherwise.
"""
function domain_is_periodic(p::ParamsBiotSavart)
    Ls = periods(p)
    !is_open_domain(Ls)
end

_extra_params(α::Zero; Ns = nothing, rcut = ∞) = (; Ns = (0, 0, 0), rcut,)  # Ns is always (0, 0, 0), no matter the input
_extra_params(α::Real; Ns, rcut = 4 / α) = (; Ns, rcut,)  # Ns is required in this case

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

function Base.summary(io::IO, p::ParamsBiotSavart{T}) where {T}
    # Print a few physical parameters (and α)
    (; Γ, a, Δ, α,) = p.common
    print(io, "ParamsBiotSavart{$T}(Γ = $Γ, a = $a, Δ = $Δ, α = $α, …)")
end

# Returns the expected period of a small-amplitude Kelvin wave of wavelength λ.
# Can be useful when setting a simulation timestep.
kelvin_wave_period(λ::Real; a, Δ, Γ) = 2 * λ^2 / Γ / (
    log(λ / (π * a)) + 1/2 - (Δ + MathConstants.γ)
)
kelvin_wave_period(p::ParamsBiotSavart, λ::Real) = kelvin_wave_period(λ; a = p.a, Δ = p.Δ, Γ = p.Γ)
