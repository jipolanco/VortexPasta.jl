function update_values_at_point_vortices!(
        fields::NamedTuple, vortices::PointVortexState, t::Real, iter::PointVortexSolver,
    )
    (; to) = iter
    # (0) Apply affect_t! function (optional)
    if iter.affect_t! !== default_callback()
        @timeit to "Affect! (fine-grained)" iter.affect_t!(iter, t)
    end
    # (1) Compute velocities and optionally streamfunction values
    @timeit to "Biot-Savart (full)" BiotSavart.compute_on_nodes!(fields, iter.cache_bs, vortices)
    # add_external_fields!(fields, iter, fs, t, iter.to)
    if iter.forcing !== NoForcing() && haskey(fields, :velocity)
        @timeit to "Add forcing" apply_forcing!(fields, iter, vortices, t)
    end
    fields
end

# Case where only the velocity is passed (generally used in internal RK substeps).
# Note: `points` only contains vortex coordinates, it doesn't contain their signs.
# But signs haven't changed from those in iter.vortices (at the start of the timestep).
function update_values_at_point_vortices!(
        vL::AbstractVector, points::AbstractVector, t::Real, iter::PointVortexSolver,
    )
    fields = (velocity = vL,)
    (; signs) = iter.vortices
    @assert eachindex(signs) == eachindex(points)
    vortices = PointVortexState(points, signs)
    update_values_at_point_vortices!(fields, vortices, t, iter)
    vL
end

function apply_forcing!(fields::NamedTuple, iter::PointVortexSolver, vortices::PointVortexState, time)
    @assert iter.forcing !== NoForcing()
    (; to) = iter
    (; vs, vf) = iter.quantities

    # At this point, fields.velocity is expected to have the self-induced velocity vs.
    # We copy that velocity into iter.vs before adding the forcing velocity to fields.velocity.
    @assert eachindex(vs) == eachindex(fields.velocity)
    @assert eachindex(vf) == eachindex(fields.velocity)
    @assert vs !== fields.velocity
    copyto!(vs, fields.velocity)

    Forcing.evaluate!(iter.forcing, iter.forcing_cache, vf, vortices, time, iter.cache_bs; to)
    @. fields.velocity = fields.velocity + vf

    fields
end
