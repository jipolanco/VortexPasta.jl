"""
    VortexFilamentProblem(fs::AbstractVector{<:AbstractFilament}, tsim::Real, p::ParamsBiotSavart)
    VortexFilamentProblem(fs::AbstractVector{<:AbstractFilament}, tspan::NTuple{2,Real}, p::ParamsBiotSavart)
    VortexFilamentProblem(checkpoint::LoadedCheckpoint, tsim::Real, p::ParamsBiotSavart)

Define a vortex filament problem.

There are basically two ways of defining a problem:

1. from scratch, from a new set of vortex filaments `fs`;

2. from a previous simulation, from a simulation state (`checkpoint`) returned from [`load_checkpoint`](@ref).

## Possible arguments

- `fs`: initial vortex positions;

- `checkpoint`: a simulation state returned from [`load_checkpoint`](@ref);

- `tsim`: total simulation time;

- `tspan = (t_begin, t_end)`: time span;

- `p`: Biot–Savart parameters (see [`ParamsBiotSavart`](@ref)).

See [`init`](@ref) for initialising a solver from a `VortexFilamentProblem`.
"""
struct VortexFilamentProblem{
        T <: AbstractFloat,
        Filaments <: VectorOfVectors{Vec3{T}, <:AbstractFilament{T}},
        Params <: ParamsBiotSavart{T},
        State <: NamedTuple,
    } <: AbstractProblem
    fs     :: Filaments
    tspan  :: NTuple{2, T}
    p      :: Params
    state  :: State  # optional simulation state, used in restarts

    # This variant assumes that types are correct (will error if they're not).
    function VortexFilamentProblem(::Type{T}, fs, tspan, p, state) where {T}
        new{T, typeof(fs), typeof(p), typeof(state)}(fs, tspan, p, state)
    end
end

function VortexFilamentProblem(
        fs_in::VectorOfFilaments,
        tspan_in::NTuple{2, Real},
        p_in::ParamsBiotSavart,
        state::NamedTuple = (;),  # empty state by default
    )
    fs = convert(VectorOfVectors, fs_in)
    # Convert all float types to the precision used to describe filaments.
    T = eltype(eltype(eltype(fs)))
    @assert T <: AbstractFloat
    tspan = convert(NTuple{2, T}, tspan_in)
    p = convert(T, p_in)
    VortexFilamentProblem(T, fs, tspan, p, state)
end

function VortexFilamentProblem(fs::VectorOfFilaments, tsim::Real, p::ParamsBiotSavart)
    tspan = (zero(tsim), tsim)
    VortexFilamentProblem(fs, tspan, p)
end

function Base.show(io::IO, prob::VortexFilamentProblem)
    (; fs, tspan, p,) = prob
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "VortexFilamentProblem with fields:")
    print(io, "\n$(prefix)├─ p: ")
    summary(io, p)
    print(io, "\n$(prefix)├─ tspan: ", tspan, " -- simulation timespan")
    _print_summary(io, fs; pre = "\n$(prefix)└─ fs: ", post = "vortex filaments at t = $(tspan[1])")
end

# This generates a summarised (shorter) version of typeof(fs).
function _typeof_summary(io::IO, fs::VectorOfVectors)
    print(io, "VectorOfVectors")
end

function _print_summary(io::IO, fs::VectorOfVectors; pre = nothing, post = nothing)
    pre === nothing || print(io, pre)
    print(io, length(fs), "-element ")
    _typeof_summary(io, fs)
    post === nothing || print(io, " -- ", post)
    nothing
end
