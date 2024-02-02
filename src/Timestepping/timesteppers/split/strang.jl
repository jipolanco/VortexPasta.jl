export Strang

"""
    Strang([fast = Midpoint()], [slow = Midpoint()]) <: SplittingScheme

2nd order Strang splitting scheme.

Uses one scheme for advancing the "fast" terms, and possibly a different scheme for the
"slow" terms. By default both schemes are taken to be the 2nd order [`Midpoint`](@ref)
method.

See [`SplittingScheme`](@ref) for more details.
"""
struct Strang{FastScheme <: TemporalScheme, SlowScheme <: TemporalScheme} <: SplittingScheme
    fast :: FastScheme
    slow :: SlowScheme
end

Strang(fast::TemporalScheme) = Strang(fast, Midpoint())
Strang() = Strang(Midpoint())

function Base.show(io::IO, scheme::Strang)
    print(io, nameof(typeof(scheme)), '(', scheme.fast, ", ", scheme.slow, ')')
end

nbuf_filaments(scheme::Strang) = 1 + max(nbuf_filaments(scheme.fast), nbuf_filaments(scheme.slow))
nbuf_velocities(scheme::Strang) = 1 + max(nbuf_velocities(scheme.fast), nbuf_velocities(scheme.slow))

function _update_velocities!(
        scheme::Strang, rhs_full!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = iter.fs, vs = iter.vs,
    ) where {F, G}
    (; fc, vc,) = cache

    ftmp = fc[1]
    vtmp = vc[1]

    cache_fast = let scheme = scheme.fast
        TemporalSchemeCache(
            scheme,
            ntuple(j -> fc[1 + j], Val(nbuf_filaments(scheme))),
            ntuple(j -> vc[1 + j], Val(nbuf_velocities(scheme))),
        )
    end
    cache_slow = let scheme = scheme.slow
        TemporalSchemeCache(
            scheme,
            ntuple(j -> fc[1 + j], Val(nbuf_filaments(scheme))),
            ntuple(j -> vc[1 + j], Val(nbuf_velocities(scheme))),
        )
    end

    gen_rhs(component) = (vs, fs, t, iter) -> rhs_full!(vs, fs, t, iter; component)

    copy!(ftmp, fs)        # initial condition for stage 1

    # Note: as opposed to other schemes, here it's not a good idea to reuse the computation
    # of the full velocity at the beginning of the timestep, since we start with the fast
    # dynamics. This means that there's an extra computation of the "slow" velocity which
    # could be avoided...

    # 1. Advance fast dynamics: t -> t + dt/2
    let component = Val(:fast), cache = cache_fast, t = t, dt = dt/2
        local rhs! = gen_rhs(component)
        rhs!(vtmp, ftmp, t, iter)
        update_velocities!(
            rhs!, advect!, cache, iter;
            resize_cache = false, t, dt, fs = ftmp, vs = vtmp,
        )
        advect!(ftmp, vtmp, dt; fbase = ftmp)
    end

    # 2. Advance slow dynamics: t -> t + dt
    let component = Val(:slow), cache = cache_slow, t = t, dt = dt
        local rhs! = gen_rhs(component)
        rhs!(vtmp, ftmp, t, iter)
        update_velocities!(
            rhs!, advect!, cache, iter;
            resize_cache = false, t, dt, fs = ftmp, vs = vtmp,
        )
        advect!(ftmp, vtmp, dt; fbase = ftmp)
    end

    # 3. Advance fast dynamics: t + dt/2 -> t + dt
    let component = Val(:fast), cache = cache_fast, t = t + dt/2, dt = dt/2
        local rhs! = gen_rhs(component)
        rhs!(vtmp, ftmp, t, iter)
        update_velocities!(
            rhs!, advect!, cache, iter;
            resize_cache = false, t, dt, fs = ftmp, vs = vtmp,
        )
        advect!(ftmp, vtmp, dt; fbase = ftmp)
    end

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end
