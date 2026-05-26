export Strang

"""
    Strang([fast = RK4()], [slow = Midpoint()]; nsubsteps::Int = 1) <: SplittingScheme

2nd order Strang splitting scheme.

Uses one scheme for advancing the "fast" terms (assumed to be cheap to compute as well), and
possibly a different scheme for the "slow" (and expensive) terms. By default these schemes
are respectively taken to be [`RK4`](@ref) and the 2nd order [`Midpoint`](@ref) method.

By default, according to Strang splitting, the fast term is advanced with a timestep of
`dt/2` (twice in a full timestep). One can pass `nsubsteps` to use even smaller timesteps
for the fast term.

See [`SplittingScheme`](@ref) for more details.
"""
struct Strang{FastScheme <: TemporalScheme, SlowScheme <: TemporalScheme} <: SplittingScheme
    fast :: FastScheme
    slow :: SlowScheme
    nsubsteps :: Int
end

Strang(fast::TemporalScheme, slow::TemporalScheme; nsubsteps::Int = 1) = Strang(fast, slow, nsubsteps)
Strang(fast::TemporalScheme; kwargs...) = Strang(fast, Midpoint(); kwargs...)
Strang(; kwargs...) = Strang(RK4(); kwargs...)

function Base.show(io::IO, scheme::Strang)
    print(io, nameof(typeof(scheme)), '(', scheme.fast, ", ", scheme.slow, "; nsubsteps = ", scheme.nsubsteps, ')')
end

nbuf_filaments(scheme::Strang) = 1 + max(nbuf_filaments(scheme.fast), nbuf_filaments(scheme.slow))
nbuf_velocities(scheme::Strang) = 1 + max(nbuf_velocities(scheme.fast), nbuf_velocities(scheme.slow))

function _update_velocities!(
        scheme::Strang, vs, rhs_full!::F, advect!::G, cache, iter::AbstractSolver;
        t = get_t(iter), dt = get_dt(iter), fs = iter.fs,
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

    function _advance_fast!(τ_start, c)
        local component = Val(:fast)
        local cache = cache_fast
        local rhs! = gen_rhs(component)
        local (; nsubsteps,) = scheme
        local τ = τ_start
        local dτ = c * dt / nsubsteps  # timestep in each substep
        for _ ∈ 1:nsubsteps
            rhs!(vtmp, ftmp, τ, iter)
            update_velocities!(
                vtmp, rhs!, advect!, cache, iter;
                resize_cache = false, t = τ, dt = dτ, fs = ftmp,
            )
            advect!(ftmp, vtmp, dτ; fbase = ftmp)
            τ += dτ
        end
        @assert τ ≈ τ_start + c * dt
    end

    function _advance_slow!(τ, c)
        local component = Val(:slow)
        local cache = cache_slow
        local rhs! = gen_rhs(component)
        local dτ = c * dt
        rhs!(vtmp, ftmp, τ, iter)
        update_velocities!(
            vtmp, rhs!, advect!, cache, iter;
            resize_cache = false, t = τ, dt = dτ, fs = ftmp,
        )
        advect!(ftmp, vtmp, dτ; fbase = ftmp)
    end

    # 1. Advance fast dynamics: t -> t + dt/2
    _advance_fast!(t, 1/2)

    # 2. Advance slow dynamics: t -> t + dt
    _advance_slow!(t, 1)

    # 3. Advance fast dynamics: t + dt/2 -> t + dt
    _advance_fast!(t + dt/2, 1/2)

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end
