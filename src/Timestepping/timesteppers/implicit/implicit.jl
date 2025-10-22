"""
    ImplicitScheme <: TemporalScheme

Abstract type representing an implicit (possibly multi-stage) scheme.

These kinds of schemes are mainly meant to be used as inner schemes (for resolving
fast-evolving dynamics which are also cheap to compute) when using
[`MultirateScheme`](@ref).
"""
abstract type ImplicitScheme <: TemporalScheme end

function mean_implicit_error(ftmp, fbase, vtmp, dt)
    ε = zero(number_type(ftmp))
    n = 0
    for (f, g, v) in zip(ftmp, fbase, vtmp)
        @inbounds for i in eachindex(f, g, v)
            δ = (f[i] - g[i]) - v[i] * dt
            ε += sum(abs2, δ)  # units of L²
            n += 1
        end
    end
    sqrt(ε / n) / dt  # units of velocity
end

function velocity_norm(vs)
    vnorm = zero(number_type(vs))
    n = 0
    for v in vs
        @inbounds for i in eachindex(v)
            vnorm += sum(abs2, v[i])
            n += 1
        end
    end
    sqrt(vnorm / n)
end

# Solve nonlinear system using fixed-point iterations (kind of brute-force method).
# We stop when the difference compared to the previous positions converges to a
# constant.
function solve_fixed_point!(
        ftmp, rhs!::F, advect!::G, iter::AbstractSolver, vtmp, v_explicit;
        dt = iter.time.dt, t = iter.time.t, cdt = dt, fbase, aI_diag,
        nmax = 200, rtol = 1e-6,
    ) where {F <: Function, G <: Function}
    @assert fbase !== ftmp
    # Currently ftmp contains an initial guess for the filament positions at time t + cdt.
    # Compute (implicit part of) RHS at the initial guess.
    rhs!(vtmp, ftmp, t + cdt, iter)
    vnorm = velocity_norm(vtmp)  # just to have a reference velocity for relative tolerance
    n = 0
    while n < nmax
        # Compute new full velocity estimate
        @. vtmp = v_explicit + aI_diag * vtmp  # full velocity estimate
        # Exit if error is small
        ε_vel = mean_implicit_error(ftmp, fbase, vtmp, dt)
        if ε_vel < rtol * vnorm
            break
        end
        n += 1
        # Update guess for filament location
        advect!(ftmp, vtmp, dt; fbase)
        rhs!(vtmp, ftmp, t + cdt, iter)  # compute RHS at the latest location
    end
    if n == nmax
        @warn lazy"reached maximum number of fixed-point iterations (n = $nmax). Convergence of implicit solver is not ensured."
    end
    vtmp  # total effective advecting velocity
end

include("CrankNicolson.jl")
