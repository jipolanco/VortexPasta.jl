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
## (2) Apply mutual friction-like forcing (normal fluid)

function apply_forcing!(fields::NamedTuple, iter::VortexFilamentSolver, fs, time, to)
    if haskey(fields, :velocity)
        # At this point, fields.velocity is expected to have the self-induced velocity vs.
        # After applying a normal fluid forcing, fields.velocity will contain the actual
        # vortex velocity vL. The self-induced velocity will be copied to iter.quantities.vs.
        _apply_forcing!(fields.velocity, iter.forcing, iter.forcing_cache, iter, fs, time, to)
    end
    fields
end

_apply_forcing!(vL, forcing::Nothing, cache::Nothing, args...) = nothing  # do nothing

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
