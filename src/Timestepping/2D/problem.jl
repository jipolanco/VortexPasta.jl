"""
    PointVortexProblem{T <: AbstractFloat}

Contains the definition of a point vortex problem.

# Construction

    PointVortexProblem(vortices::PointVortexState, tspan::NTuple{2, Real}, p::ParamsPointVortex)
    PointVortexProblem(vortices::PointVortexState, tsim::Real, p::ParamsPointVortex)

## Possible arguments

- `vortices`: initial vortex positions and signs (see [`PointVortexState`](@ref));

- `tsim`: total simulation time;

- `tspan = (t_begin, t_end)`: time span;

- `p`: point vortex parameters (see [`ParamsPointVortex`](@ref)).

Once constructed, the problem can be passed to [`init`](@ref) to set-up a simulation.

See also [`VortexFilamentProblem`](@ref) for the 3D version of this function.
"""
struct PointVortexProblem{
        T <: AbstractFloat,
        Vortices <: PointVortexState{T},
        Params <: ParamsPointVortex{T},
        State <: NamedTuple,
    } <: AbstractProblem
    vortices::Vortices
    tspan::NTuple{2, T}
    p::Params
    state::State  # optional simulation state, used in restarts

    # This variant assumes that types are correct (will error if they're not).
    function PointVortexProblem(::Type{T}, vortices::PointVortexState{T}, tspan, p::ParamsPointVortex{T}, state) where {T}
        return new{T, typeof(vortices), typeof(p), typeof(state)}(vortices, tspan, p, state)
    end
end

function PointVortexProblem(
        vortices::PointVortexState{T},
        tspan_in::NTuple{2, Real},
        p_in::ParamsPointVortex,
        state::NamedTuple = (;),  # empty state by default
    ) where {T <: AbstractFloat}
    tspan = convert(NTuple{2, T}, tspan_in)
    p = convert(T, p_in)
    PointVortexProblem(T, vortices, tspan, p, state)
end

function PointVortexProblem(vortices::PointVortexState, tsim::Real, p::ParamsPointVortex)
    tspan = (zero(tsim), tsim)
    PointVortexProblem(vortices, tspan, p)
end

function Base.show(io::IO, prob::PointVortexProblem)
    (; vortices, tspan, p,) = prob
    prefix = get(io, :prefix, " ")  # single space by default
    print(io, "PointVortexProblem with fields:")
    print(io, "\n$(prefix)├─ p: ")
    summary(io, p)
    print(io, "\n$(prefix)├─ tspan: ", tspan, " -- simulation timespan")
    print(io, "\n$(prefix)└─ vortices: ", vortices, " at t = ", tspan[1])
    nothing
end
