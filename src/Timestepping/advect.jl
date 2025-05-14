# These functions are used to advect filament locations in time: t -> t + δt, according to a
# computed (mean) velocity throughout the timestep.
#
# There are two possible cases:
#
# (1) δt is the full simulation timestep Δt, and we're computing the positions at time
#     t + Δt from positions at time t and the mean velocities obtained from the chosen
#     temporal scheme;
#
# (2) δt is a fraction of Δt, e.g. in the case of a RK substep.

# This function will dispatch to one of the two possible cases, according to the value of fbase:
#
# - `nothing` (default), used when the filaments `fs` should be advected with the
#   chosen velocity at the end of a timestep (case 1);
#
# - some list of filaments similar to `fs`, used to advect from `fbase` to `fs` with the
#   chosen velocity (case 2). This is generally done from within RK stages, and in this case
#   we do not recompute the filament parametrisation (knots).
advect_filaments!(fs, vL, dt; fbase = nothing, kws...) =
    _advect_filaments!(fs, fbase, vL, dt; kws...)

# Variant called at the end of a timestep (case 1)
function _advect_filaments!(fs, fbase::Nothing, vL, dt)
    for i ∈ eachindex(fs, vL)
        @inbounds _advect_filament!(fs[i], nothing, vL[i], dt)
    end
    fs
end

function _advect_filament!(
        f::AbstractFilament, fbase::Nothing, vL::VectorOfVelocities, dt::Real,
    )
    Xs = nodes(f)
    for i ∈ eachindex(Xs, vL)
        @inbounds Xs[i] = Xs[i] + dt * vL[i]
    end
    Filaments.update_coefficients!(f)  # note: in this case we recompute the knots
    nothing
end

# Variant called from within RK stages (case 2)
function _advect_filaments!(fs::T, fbase::T, vL, dt) where {T}
    for i ∈ eachindex(fs, fbase, vL)
        @inbounds _advect_filament!(fs[i], fbase[i], vL[i], dt)
    end
    fs
end

function _advect_filament!(f::T, fbase::T, vL, dt) where {T <: AbstractFilament}
    Xs = nodes(f)
    Ys = nodes(fbase)
    for i ∈ eachindex(Xs, Ys, vL)
        @inbounds Xs[i] = Ys[i] + dt * vL[i]
    end
    # Compute interpolation coefficients making sure that knots are preserved (and not
    # recomputed from new locations).
    Filaments.update_coefficients!(f; knots = knots(fbase))
    nothing
end
