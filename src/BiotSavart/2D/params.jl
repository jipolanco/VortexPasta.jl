struct ParamsCommon2D{T <: AbstractFloat, Splitting <: PeriodicSplitting{2}} <: AbstractParamsCommon{T}
    Γ::T                  # vortex circulation
    splitting::Splitting  # Ewald splitting kernel (Gaussian, ...)
    Ls::NTuple{2, T}      # size of unit cell (= period in each direction)
    function ParamsCommon2D{T}(Γ, splitting_in, Ls_in) where {T}
        splitting = convert_floats(T, splitting_in)
        Ls = convert.(T, Ls_in)
        new{T, typeof(splitting)}(Γ, splitting, Ls)
    end
end

function Base.convert(::Type{T}, p::ParamsCommon2D) where {T <: AbstractFloat}
    (; Γ, splitting, Ls) = p
    ParamsCommon2D{T}(Γ, splitting, Ls)  # converts all floats to type T
end

function Base.show(io::IO, p::ParamsCommon2D)
    (; Γ,) = p
    indent = get(io, :indent, 0)
    io_sub = IOContext(io, :indent => indent + 5)
    print(io, "\n - Vortex circulation: Γ  = ", Γ)
    print(io, "\n - Ewald splitting kernel:\n")
    print(io_sub, p.splitting)
    nothing
end

function to_hdf5(g, p::ParamsCommon2D{T}) where {T}
    g["Gamma"] = p.Γ
    g["Ls"] = collect(T, p.Ls)    # this allows converting Infinity() to a float (e.g. Inf)
    g["splitting"] = string(p.splitting)
    nothing
end

## ================================================================================ ##

struct ParamsShortRange2D{
        T <: AbstractFloat,
        Backend <: ShortRangeBackend,
        Common <: ParamsCommon2D{T},
    } <: AbstractParamsShortRange{T}
    backend::Backend
    common::Common    # common parameters (Γ, splitting, Ls)
    rcut::T           # cutoff distance
    # use_simd::Bool  # not yet supported (TODO?)
end

function Base.show(io::IO, p::ParamsShortRange2D)
    print(io, "\n   * Short-range backend:  ", p.backend)
    # print(io, "\n   * Short-range uses SIMD:  ", p.use_simd)
    nothing
end

function to_hdf5(g, p::ParamsShortRange2D)
    g["backend_shortrange"] = string(p.backend)
    # g["use_simd"] = p.use_simd
    nothing
end

# Create a new ParamsShortRange2D based on the type of ParamsCommon2D.
function change_float_type(p::ParamsShortRange2D, common::ParamsCommon2D{T}) where {T}
    ParamsShortRange2D(p.backend, common, convert(T, p.rcut))
end

backend(p::ParamsShortRange2D) = p.backend

## ================================================================================ ##

struct ParamsLongRange2D{
        T <: AbstractFloat,
        Backend <: LongRangeBackend,
        Common <: ParamsCommon2D{T},
    } <: AbstractParamsLongRange{T}
    backend::Backend
    common::Common  # common parameters (Γ, splitting, Ls)
    Ns::Dims{2}     # grid dimensions for FFTs
    truncate_spherical::Bool  # if true, perform spherical truncation in Fourier space
end

has_real_to_complex(p::ParamsLongRange2D) = has_real_to_complex(p.backend)
maximum_wavenumber(p::ParamsLongRange2D) = maximum_wavenumber(p.Ns, p.common.Ls)

function Base.show(io::IO, p::ParamsLongRange2D{T}) where {T}
    print(io, "\n   * Long-range backend:   ", p.backend)
    print(io, "\n   * Long-range spherical truncation: ", p.truncate_spherical)
    nothing
end

function to_hdf5(g, p::ParamsLongRange2D{T}) where {T}
    (; Ns,) = p
    g["backend_longrange"] = string(p.backend)
    g["Ns"] = collect(Ns)
    g["truncate_spherical"] = p.truncate_spherical
    nothing
end

function change_float_type(p::ParamsLongRange2D, common::ParamsCommon2D{T}) where {T}
    ParamsLongRange2D(p.backend, common, p.Ns, p.truncate_spherical)
end

backend(p::ParamsLongRange2D) = p.backend

## ================================================================================ ##

"""
    ParamsPointVortex{T <: AbstractFloat}

Contains parameters for point vortex computations in periodic systems using fast Ewald splitting.

This is the 2D version of [`ParamsBiotSavart`](@ref).

The type parameter `T` corresponds to the precision used in computations
(typically `Float64` or `Float32`).

# Construction

    ParamsPointVortex([T = Float64]; Γ, splitting, optional_kws...)

where the optional parameter `T` sets the numerical precision.

Mandatory and optional keyword arguments are detailed in the extended help below.

# Extended help

## Mandatory keyword arguments

- `Γ::Real`: vortex circulation (assumed constant);

- `splitting::PeriodicSplitting`: a choice of Ewald splitting kernel for FFT-based
  long-range computations.
  Possible choices are:
  - [`GaussianSplitting`](@ref): standard Gaussian splitting window;
  - [`KaiserBesselSplitting`](@ref): a more efficient splitting window (**not yet implemented in 2D**).
  Note that it is currently not possible to disable periodicity in 2D point vortex computations.

## Optional keyword arguments and their defaults

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

"""
struct ParamsPointVortex{
        T <: AbstractFloat,
        Common <: ParamsCommon2D{T},
        ShortRange <: ParamsShortRange2D,
        LongRange <: ParamsLongRange2D,
    } <: AbstractParams{T}
    common     :: Common
    shortrange :: ShortRange
    longrange  :: LongRange
end

function ParamsPointVortex(
        Γ::T, splitting::PeriodicSplitting{2};
        backend_short::ShortRangeBackend = default_short_range_backend(splitting),
        backend_long::LongRangeBackend = default_long_range_backend(splitting),
        longrange_truncate_spherical::Bool = false,
    ) where {T <: AbstractFloat}
    Ls = periods(splitting)::NTuple{2, Real}
    Ns = fourier_grid_size(splitting)
    rcut = T(cutoff_distance(splitting))
    common = ParamsCommon2D{T}(Γ, splitting, Ls)
    sr = ParamsShortRange2D(backend_short, common, rcut)
    lr = ParamsLongRange2D(backend_long, common, Ns, longrange_truncate_spherical)
    ParamsPointVortex(common, sr, lr)
end

ParamsPointVortex(Γ::T; splitting::PeriodicSplitting{2}, kws...) where {T <: AbstractFloat} =
    ParamsPointVortex(Γ, splitting; kws...)
ParamsPointVortex(Γ::Real; kws...) = ParamsPointVortex(Float64(Γ); kws...)  # e.g. if we pass Γ = 1 (integer)
ParamsPointVortex(; Γ, kws...) = ParamsPointVortex(Γ; kws...)

function Base.convert(::Type{T}, p::ParamsPointVortex) where {T <: AbstractFloat}
    common = convert(T, p.common)
    shortrange = change_float_type(p.shortrange, common)
    longrange = change_float_type(p.longrange, common)
    ParamsPointVortex(common, shortrange, longrange)
end

Base.ndims(::ParamsPointVortex) = 2

circulation(p::ParamsPointVortex) = p.common.Γ
periods(p::ParamsPointVortex) = p.common.Ls
domain_is_periodic(p::ParamsPointVortex) = true  # we're always periodic in 2D

# This is for convenience: doing p.splitting is equivalent to p.common.splitting.
@inline function Base.getproperty(p::ParamsPointVortex, name::Symbol)
    common = getfield(p, :common)
    if hasproperty(common, name)
        getproperty(common, name)
    else
        getfield(p, name)
    end
end

function Base.propertynames(p::ParamsPointVortex, private::Bool = false)
    (fieldnames(typeof(p))..., propertynames(p.common, private)...)
end

function Base.show(io::IO, p::ParamsPointVortex{T}) where {T}
    print(io, "ParamsPointVortex{$T} with:")
    show(io, p.common)
    show(io, p.shortrange)
    show(io, p.longrange)
    nothing
end

function Base.summary(io::IO, p::ParamsPointVortex{T}) where {T}
    # Print a few physical and numerical parameters
    (; splitting,) = p.common
    Γ = round(p.common.Γ; sigdigits = 3)
    print(io, "ParamsPointVortex{$T}(Γ ≈ $Γ, splitting = ")
    summary(io, splitting)
    print(io, ", …)")
end
