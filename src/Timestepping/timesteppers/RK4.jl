struct RK4 <: ExplicitTemporalScheme end

struct RK4Cache{
        Velocities <: VectorOfVelocities,
    } <: TemporalSchemeCache
    vs :: NTuple{4, Velocities}
end