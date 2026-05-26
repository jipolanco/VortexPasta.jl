export Strang4

"""
    Strang([fast = RK4()], [slow = RK4()]; nsubsteps::Int = 1) <: SplittingScheme

4th order Strang splitting scheme.

Uses one scheme for advancing the "fast" terms (assumed to be cheap to compute as well), and
possibly a different scheme for the "slow" (and expensive) terms. By default, to preserve
the global 4th-order accuracy, both schemes are taken to be [`RK4`](@ref) method.

According to the Triple jump composition procedure, the time step is split into three substeps 
using coefficients (γ, 1 - 2γ, γ) with 'γ = \frac{1}{2 - 2^\frac{1}{3}}'. One can pass 
`nsubsteps` to use even smaller timesteps for the fast term.

See [`SplittingScheme`](@ref) for more details.
"""
struct Strang4{FastScheme <: TemporalScheme, SlowScheme <: TemporalScheme} <: SplittingScheme
    fast :: FastScheme
    slow :: SlowScheme
    nsubsteps :: Int
end

Strang4(fast::TemporalScheme, slow::TemporalScheme; nsubsteps::Int = 1) = Strang4(fast, slow, nsubsteps)
Strang4(fast::TemporalScheme; kwargs...) = Strang4(fast, RK4(); kwargs...)
Strang4(; kwargs...) = Strang4(RK4(); kwargs...)

function Base.show(io::IO, scheme::Strang4)
    print(io, nameof(typeof(scheme)), '(', scheme.fast, ", ", scheme.slow, "; nsubsteps = ", scheme.nsubsteps, ')')
end

nbuf_filaments(scheme::Strang4) = 1 + max(nbuf_filaments(scheme.fast), nbuf_filaments(scheme.slow))
nbuf_velocities(scheme::Strang4) = 1 + max(nbuf_velocities(scheme.fast), nbuf_velocities(scheme.slow))

function _update_velocities!(
        scheme::Strang4, vs, rhs_full!::F, advect!::G, cache, iter::AbstractSolver;
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

    γ1 = 1 / (2 - cbrt(2))
    # 1. First step: t -> t + γ * dt
    _advance_fast!(t, γ1 / 2)
    _advance_slow!(t, γ1)

    # 2. Second step: t + γ * dt -> t + (1 - 2γ) * dt
    _advance_fast!(t + γ1 / 2 * dt, (1 - γ1) / 2)
    _advance_slow!(t + γ1 * dt, (1 - 2 * γ1))

    # 3. Third step: t + (1 - 2γ) * dt -> t + γ * dt
    _advance_fast!(t + dt / 2, (1 - γ1) / 2)
    _advance_slow!(t + (1 - γ1) * dt, γ1)
    _advance_fast!(t + (2 - γ1) / 2 * dt, γ1 / 2)

    # Now ftmp is at the final position. We compute the effective velocity to go from fs to
    # ftmp (for consistency with other schemes).
    for i ∈ eachindex(fs, ftmp, vs)
        @. vs[i] = (ftmp[i] - fs[i]) / dt
    end

    vs
end