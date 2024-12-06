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
## (2) Apply mutual friction forcing (normal fluid)

function apply_forcing!(fields::NamedTuple, iter::VortexFilamentSolver, fs, time, to)
    if haskey(fields, :velocity)
        _apply_forcing!(fields.velocity, iter.forcing, iter, fs, time, to)
    end
    fields
end

_apply_forcing!(vs_all, ::Nothing, args...) = nothing  # do nothing

function _apply_forcing!(vs_all, forcing::NormalFluidForcing, iter, fs, t, to)
    vn_all = iter.vn
    tangents_all = iter.tangents
    @assert eachindex(vn_all) == eachindex(tangents_all) == eachindex(vs_all) == eachindex(fs)
    @timeit to "Add forcing" begin
        @inbounds Threads.@threads for n ∈ eachindex(fs, vs_all, vn_all, tangents_all)
            f = fs[n]
            vs = vs_all[n]
            vn = vn_all[n]
            tangents = tangents_all[n]
            # Compute tangents (TODO: can we reuse them?)
            for i ∈ eachindex(f, tangents)
                tangents[i] = f[i, UnitTangent()]
            end
            Forcing.get_velocities!(forcing, vn, f)    # evaluate normal fluid velocities
            Forcing.apply!(forcing, vs, vn, tangents)  # modify vs according to the given forcing
        end
    end
    nothing
end