"""
    LongRangeCacheState

Describes the current state of the long-range fields in Fourier space.

It has two fields, which allow to know which field is currently stored and its smoothing scale
(width ``σ`` of Gaussian filter):

- `quantity::Symbol` (can be `:undef`, `:vorticity`, `:velocity`, `:streamfunction`)

- `smoothing_scale::Float64`: width ``σ`` of applied Gaussian filter. This is typically
  0 (if the field has not been smoothed), but that can change if one calls
  [`to_coarse_grained_vorticity!`](@ref).

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
        ParamsAll <: ParamsBiotSavart{T},
        WaveNumbers <: NTuple{3, AbstractVector{T}},
        PointCharges <: PointData{T},
        OutputVectors <: NamedTuple,
        FourierVectorField <: StructArray{Vec3{Complex{T}}, 3},
        RealScalarField <: AbstractArray{T, 3},
        Timer <: TimerOutput,
    }
    params      :: Params
    params_all  :: ParamsAll  # note: params === params_all.longrange
    wavenumbers :: WaveNumbers
    pointdata   :: PointCharges        # non-uniform data in physical space
    outputs     :: OutputVectors       # output velocity and streamfunction fields (as linear vectors)
    uhat        :: FourierVectorField  # uniform Fourier-space data (3 × [Nx, Ny, Nz])
    vorticity_prefactor :: T             # prefactor Γ/V appearing in the Fourier transform of the vorticity
    ewald_gaussian    :: RealScalarField   # real-valued Ewald Gaussian smoothing in Fourier space ([Nx, Ny, Nz])
    state           :: LongRangeCacheState
    to              :: Timer             # timer for operations run on the device
end

get_wavenumbers(c::LongRangeCacheCommon) = c.wavenumbers
get_wavenumbers(c::LongRangeCache) = get_wavenumbers(c.common)

default_interpolation_output(c::LongRangeCacheCommon) = c.outputs.default
default_interpolation_output(c::LongRangeCache) = default_interpolation_output(c.common)

# Callback to be used right before interpolation into physical space, which applies the
# Gaussian filter in Ewald's method.
# We wrap it in a struct to make sure it works correctly on GPUs.
struct EwaldInterpolationCallback{ScalarField <: AbstractArray} <: Function
    data :: ScalarField
end

# Note: to be compatible with NonuniformFFTs, the callback must have the signature f(û::NTuple{3}, idx::NTuple{3}).
# (Same signature as default_callback_interp.)
@inline function (f::EwaldInterpolationCallback)(û::NTuple{3,T}, idx::NTuple{3}) where {T}
    @inbounds op = f.data[idx...]
    map(v -> T(v * op), û)
end

@inline get_ewald_interpolation_callback(c::LongRangeCacheCommon) = EwaldInterpolationCallback(c.ewald_gaussian)
@inline get_ewald_interpolation_callback(c::LongRangeCache) = get_ewald_interpolation_callback(c.common)

function LongRangeCacheCommon(
        params_all::ParamsBiotSavart{T},
        wavenumbers::NTuple{3, AbstractVector},
        pointdata_in::PointData{T},
    ) where {T}
    pcommon = params_all.common
    params = params_all.longrange
    (; Γ, Ls, α,) = pcommon
    (; backend,) = params
    @assert T === eltype(pcommon)
    @assert α !== Zero()
    Nks = map(length, wavenumbers)
    uhat = init_fourier_vector_field(backend, T, Nks) :: StructArray{Vec3{Complex{T}}, 3}
    vorticity_prefactor = Γ / prod(Ls)
    ka_backend = KA.get_backend(backend)  # CPU, CUDABackend, ROCBackend, ...
    wavenumbers = adapt(ka_backend, wavenumbers)  # copy wavenumbers onto device if needed
    pointdata = adapt(ka_backend, pointdata_in)      # create PointData replica on the device if needed
    if pointdata === pointdata_in       # basically if ka_backend isa CPU
        pointdata = copy(pointdata_in)  # make sure pointdata and pointdata_in are not aliased!
    end
    outputs = (;
        velocity = similar(pointdata.charges),
        streamfunction = similar(pointdata.charges),
        default = similar(pointdata.charges),  # this is the "default" output when no output has been selected
    )
    ewald_gaussian = init_ewald_gaussian_operator(T, backend, wavenumbers, α)
    state = LongRangeCacheState()
    to = TimerOutput()
    LongRangeCacheCommon(params, params_all, wavenumbers, pointdata, outputs, uhat, vorticity_prefactor, ewald_gaussian, state, to)
end

has_real_to_complex(c::LongRangeCacheCommon) = has_real_to_complex(c.params)
