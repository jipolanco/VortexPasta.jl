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
    if haskey(fields, :velocity)
        if fields.velocity !== iter.vs
            # Make sure iter.vs contains the newly computed velocity (may be used in forcing or dissipation)
            @timeit to "Copy self-induced velocity" copyto!(iter.vs, fields.velocity)
        end
        if iter.forcing !== NoForcing()
            @timeit to "Add forcing" apply_forcing!(fields, iter, vortices, t)
        end
        if iter.dissipation !== NoDissipation()
            @timeit to "Add dissipation" apply_dissipation!(fields, iter, vortices, t)
        end
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
    (; vf) = iter.quantities
    Forcing.evaluate!(iter.forcing, iter.forcing_cache, vf, vortices, time, iter.cache_bs; to)
    @. fields.velocity = fields.velocity + vf
    fields
end

function apply_dissipation!(fields::NamedTuple, iter::PointVortexSolver, vortices::PointVortexState, time)
    @assert iter.dissipation !== NoDissipation()
    (; to) = iter
    (; vs, vdiss) = iter.quantities

    # At this point, iter.vs is expected to already have the self-induced velocity vs.
    # If forcing is enabled, we currently have fields.velocity = iter.vs + iter.vf.
    # Moreover, at the beginning of a full timestep, fields.velocity is aliased to iter.vL.
    @assert eachindex(vs) == eachindex(fields.velocity)
    @assert eachindex(vdiss) == eachindex(fields.velocity)

    Forcing.evaluate!(
        iter.dissipation, iter.dissipation_cache, vdiss, vs, vortices, time,
        iter.cache_bs; to
    )
    @. fields.velocity = fields.velocity + vdiss

    fields
end
