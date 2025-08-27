## ================================================================================ ##
## (1) Apply external fields and stretching velocity

function add_external_fields!(fields::NamedTuple, iter::VortexFilamentSolver, fs, t, to)
    (; external_fields, stretching_velocity,) = iter
    if haskey(fields, :velocity)
        _add_external_field!(fields.velocity, external_fields.velocity, fs, t, to)
        _add_stretching_velocity!(fields.velocity, stretching_velocity, fs, to)
    end
    if haskey(fields, :streamfunction)
        # We multiply the external streamfunction by 2 to get the right kinetic energy.
        _add_external_field!(fields.streamfunction, external_fields.streamfunction, fs, t, to; factor = 2)
    end
    fields
end

_add_external_field!(vs_all, ::Nothing, args...; kws...) = vs_all  # do nothing

function _add_external_field!(vs_all, vext::F, fs, time, to; factor = 1) where {F <: Function}
    @assert eachindex(vs_all) == eachindex(fs)
    @timeit to "Add external field" begin
        for (f, vs) ∈ zip(fs, vs_all)
            for i ∈ eachindex(vs, f)
                @inbounds vs[i] = vs[i] + factor * vext(f[i], time)
            end
        end
    end
    vs_all
end

_add_stretching_velocity!(vs_all, ::Nothing, args...) = vs_all  # do nothing

function _add_stretching_velocity!(vs_all, stretching_velocity::F, fs, to) where {F <: Function}
    @assert eachindex(vs_all) == eachindex(fs)
    @timeit to "Add stretching velocity" begin
        for (f, vs) ∈ zip(fs, vs_all)
            @inbounds for i ∈ eachindex(vs, f)
                ρ⃗ = f[i, CurvatureVector()]
                ρ = sqrt(sum(abs2, ρ⃗))  # = |ρ⃗|
                n̂ = ρ⃗ ./ ρ  # normal vector
                vs[i] = vs[i] - n̂ * stretching_velocity(ρ)
            end
        end
    end
    vs_all
end

## ================================================================================ ##
## (2) Apply optional external forcing and dissipation terms.

function apply_forcing!(fields::NamedTuple, iter::VortexFilamentSolver, fs, time, to)
    if haskey(fields, :velocity)
        # At this point, fields.velocity is expected to have the self-induced velocity vs.
        # After applying a normal fluid forcing, fields.velocity will contain the actual
        # vortex velocity vL. The self-induced velocity will be copied to iter.quantities.vs.
        _apply_forcing!(fields.velocity, iter.forcing, iter.forcing_cache, iter, fs, time, to)
        # NOTE: dissipation should be applied _after_ forcing, since forcing assumes that vL
        # currently contains the Biot-Savart velocity (with no extra terms).
        _apply_dissipation!(fields.velocity, iter.dissipation, iter.dissipation_cache, iter, fs, time, to)
    end
    fields
end

# For consistency with other forcing methods, we must copy vL -> iter.vs (the self-induced BS velocity).
# This is needed in particular if we're including a dissipation term.
function _apply_forcing!(vL, forcing::NoForcing, cache, iter, args...)
    # The condition below basically means that there is a dissipation term.
    # Note: iter.vL is not necessarily the same as vL, especially if we're inside of a RK substep.
    if iter.vs !== iter.vL
        copyto!(iter.vs, vL)
    end
    nothing
end

_apply_dissipation!(vL, dissipation::NoDissipation, cache, args...) = nothing  # do nothing

# Note: the cache is currently not used by NormalFluidForcing (it's empty anyway)
function _apply_forcing!(vL_all, forcing::NormalFluidForcing, cache, iter, fs, t, to)
    # Note: inside a RK substep, quantities.{vs,vn,tangents} are used as temporary buffers
    # (and in fact we don't really need vs).
    (; quantities,) = iter
    vs_all = quantities.vs  # self-induced velocities will be copied here
    vn_all = quantities.vn  # normal fluid velocities will be computed here
    tangents_all = quantities.tangents  # local tangents (already computed)
    @assert vs_all !== vL_all
    @assert eachindex(vL_all) == eachindex(vn_all) == eachindex(tangents_all) == eachindex(vs_all) == eachindex(fs)
    scheduler = DynamicScheduler()  # for threading
    @timeit to "Add forcing" begin
        @inbounds for n ∈ eachindex(fs)
            f = fs[n]
            vs = vs_all[n]
            vL = vL_all[n]
            vn = vn_all[n]
            tangents = tangents_all[n]
            # At input, vL contains the self-induced velocity vs.
            # We copy vL -> vs before modifying vL with the actual vortex velocities.
            tforeach(eachindex(f, tangents); scheduler) do i
                vs[i] = vL[i]  # usually this is the Biot-Savart velocity
                vn[i] = forcing.vn(f[i])  # evaluate normal fluid velocity
            end
            Forcing.apply!(forcing, vL, vn, tangents; scheduler)  # compute vL according to the given forcing (vL = vs at input)
        end
    end
    nothing
end

function _apply_forcing!(vL_all, forcing::FourierBandForcing, cache, iter, fs, t, to)
    @assert eachindex(vL_all) === eachindex(fs)
    (; quantities,) = iter
    vs_all = quantities.vs  # self-induced velocities will be copied here
    v_ns_all = quantities.v_ns  # slip velocity used in forcing
    Forcing.update_cache!(cache, forcing, iter.cache_bs)
    scheduler = DynamicScheduler()  # for threading
    @timeit to "Add forcing" begin
        @inbounds for n in eachindex(fs)
            f = fs[n]
            vL = vL_all[n]  # currently contains self-induced velocity vs
            copyto!(vs_all[n], vL)
            # Compute vL according to the given forcing (vL = vs at input).
            # This will also store filtered slip velocities in v_ns.
            Forcing.apply!(forcing, cache, vL, f; vdiff = v_ns_all[n], scheduler)
        end
    end
    nothing
end

function _apply_forcing!(vL_all, forcing::FourierBandForcingBS, cache, iter, fs, t, to)
    @assert eachindex(vL_all) === eachindex(fs)
    (; quantities,) = iter
    vs_all = quantities.vs  # self-induced velocities will be copied here
    vf_all = quantities.vf  # forcing velocities will be copied here
    tangents_all = quantities.tangents  # local tangents (already computed)
    Forcing.update_cache!(cache, forcing, iter.cache_bs)
    scheduler = DynamicScheduler()  # for threading
    @timeit to "Add forcing" begin
        ε_total = zero(number_type(vs_all))
        @inbounds for n in eachindex(fs)
            f = fs[n]
            vL = vL_all[n]  # currently contains self-induced velocity vs
            copyto!(vs_all[n], vL)
            # Compute vL according to the given forcing (vL = vs at input).
            ε_total += Forcing.apply!(forcing, cache, vL, f; scheduler)
        end
        if forcing.ε_target != 0 && ε_total != 0
            # In this case, vL only contains the forcing velocities, which need to be
            # rescaled to get the wanted energy injection rate (estimated).
            @assert forcing.α == 0
            α = forcing.ε_target / ε_total
            @inbounds for n in eachindex(fs, vs_all, vL_all)
                f = fs[n]
                vL = vL_all[n]
                vs = vs_all[n]
                tangents = tangents_all[n]
                if forcing.α′ == 0
                    for i in eachindex(f, vs, vL)
                        vL[i] = vs[i] + α * vL[i]
                    end
                else
                    for i in eachindex(f, vs, vL, tangents)
                        s⃗′ = tangents[i]
                        vL[i] = vs[i] + α * vL[i] - forcing.α′ * (s⃗′ × vL[i])
                    end
                end
            end
        end
        @. vf_all = vL_all - vs_all
    end
    nothing
end

_apply_dissipation!(vL_all, dissipation::NoDissipation, cache, iter, fs, t, to) = nothing

function _apply_dissipation!(vL_all, dissipation::AbstractDissipation, cache, iter, fs, t, to)
    @assert eachindex(vL_all) === eachindex(fs)
    (; quantities,) = iter
    @assert quantities.vs !== vL_all  # not aliased
    scheduler = DynamicScheduler()  # for threading
    @timeit to "Add dissipation" begin
        # This will:
        # - add dissipation term to vL_all
        # - write dissipation term into vdiss
        Forcing.apply!(
            dissipation, cache,
            vL_all, quantities.vdiss, quantities.vs, fs;
            scheduler,
        )
    end
    nothing
end
