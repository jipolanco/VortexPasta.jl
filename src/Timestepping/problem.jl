"""
    VortexFilamentProblem
    VortexFilamentProblem(fs::AbstractVector{<:AbstractFilament}, tspan::NTuple{2}, p::ParamsBiotSavart)

Define a vortex filament problem.

Arguments:

- `fs`: initial vortex positions;

- `tspan = (t_begin, t_end)`: time span;

- `p`: Biot–Savart parameters (see [`ParamsBiotSavart`](@ref)).

See [`init`](@ref) for initialising a solver from a `VortexFilamentProblem`.
"""
struct VortexFilamentProblem{
        T <: AbstractFloat,
        Filaments <: VectorOfVectors{Vec3{T}, <:AbstractFilament{T}},
        Params <: ParamsBiotSavart{T},
    } <: AbstractProblem
    fs    :: Filaments
    tspan :: NTuple{2, T}
    p     :: Params

    # This variant assumes that types are correct (will error if they're not).
    function VortexFilamentProblem(::Type{T}, fs, tspan, p) where {T}
        new{T, typeof(fs), typeof(p)}(fs, tspan, p)
    end
end

function VortexFilamentProblem(
        fs_in::VectorOfFilaments,
        tspan_in::NTuple{2, Real},
        p_in::ParamsBiotSavart,
    )
    fs = convert(VectorOfVectors, fs_in)
    # Convert all float types to the precision used to describe filaments.
    T = eltype(eltype(eltype(fs)))
    @assert T <: AbstractFloat
    tspan = convert(NTuple{2, T}, tspan_in)
    p = convert(T, p_in)
    VortexFilamentProblem(T, fs, tspan, p)
end

function Base.show(io::IO, prob::VortexFilamentProblem)
    (; fs, tspan, p,) = prob
    print(io, "VortexFilamentProblem with fields:")
    print(io, "\n - `p`: ")
    summary(io, p)
    print(io, "\n - `tspan`: ", tspan)
    _print_summary(io, fs; pre = "\n - `fs`: ")
end

# This generates a summarised (shorter) version of typeof(fs).
function _typeof_summary(io::IO, fs::VectorOfVectors)
    V = eltype(fs)  # e.g. ClosedFilament{Float64, ...}
    T = eltype(V)   # e.g. Vec3{Float64}
    print(io, "VectorOfVectors{$T}")
end

function _print_summary(io::IO, fs::VectorOfVectors; pre = nothing, post = nothing)
    pre === nothing || print(io, pre)
    print(io, length(fs), "-element ")
    _typeof_summary(io, fs)
    post === nothing || print(io, post)
    nothing
end
