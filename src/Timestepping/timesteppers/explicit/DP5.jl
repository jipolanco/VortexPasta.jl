export DP5

"""
    DP5 <: ExplicitScheme

Dormand–Prince embedded 4/5 Runge–Kutta method.

For now, the implementation is **experimental** and can be further optimised.

See [Wikipedia](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method) for details.
"""
struct DP5 <: ExplicitScheme end

# TODO
# - can we reduce the number of buffers?
# - take advantage of the FSAL (First Same As Last) property (see wikipedia)
# - enable error-based adaptivity for embedded methods

# Number of buffers needed to hold "intermediate" filaments and velocities.
nbuf_filaments(::DP5) = 1
nbuf_velocities(::DP5) = 7

function _update_velocities!(
        ::DP5, rhs!::F, advect!::G, cache, iter::AbstractSolver,
    ) where {F <: Function, G <: Function}
    (; fs, vs, to,) = iter
    (; fc, vc,) = cache

    t = get_t(iter)
    dt = get_dt(iter)

    fbase = fs
    ftmp = fc[1]

    cs = (0, 1/5, 3/10, 4/5, 8/9, 1, 1)

    # Stage 1
    @. vc[1] = vs

    # Stage 2
    @. vs = vc[1]  # advection velocity
    advect!(ftmp, vs, cs[2] * dt; fbase)
    rhs!(vc[2], ftmp, t + cs[2] * dt, iter)

    # Stage 3
    @. vs = (vc[1] + 3 * vc[2]) / 4
    advect!(ftmp, vs, cs[3] * dt; fbase)
    rhs!(vc[3], ftmp, t + cs[3] * dt, iter)

    # Stage 4
    @. vs = (11 * vc[1] - 42 * vc[2] + 40 * vc[3]) / 9
    advect!(ftmp, vs, cs[4] * dt; fbase)
    rhs!(vc[4], ftmp, t + cs[4] * dt, iter)

    # Stage 5
    @. vs = (
         4843 / 1458 * vc[1] +
        -3170 /  243 * vc[2] +
         8056 /  729 * vc[3] +
          -53 /  162 * vc[4]
    )
    advect!(ftmp, vs, cs[5] * dt; fbase)
    rhs!(vc[5], ftmp, t + cs[5] * dt, iter)

    # Stage 6
    @. vs = (
        9017/3168 * vc[1] +
        -355/33 * vc[2] +
        46732/5247 * vc[3] +
        49/176 * vc[4] +
        -5103/18656 * vc[5]
    )
    advect!(ftmp, vs, cs[6] * dt; fbase)
    rhs!(vc[6], ftmp, t + cs[6] * dt, iter)

    # Stage 7
    @. vs = (
        35/384 * vc[1] +
        # 0/1 * vc[2] +
        500/1113 * vc[3] +
        125/192 * vc[4] +
        -2187/6784 * vc[5] +
        11/84 * vc[6]
    )  # this velocity is supposed to be 5-th order accurate?

    # If I understand correctly, the following is only needed for error estimation:
    # advect!(ftmp, vs, cs[7] * dt; fbase)
    # rhs!(vc[7], ftmp, t + cs[7] * dt, iter)

    # Get 4-th order accurate velocity (written onto vc[1])
    # @. vc[1] = (
    #     5179/57600 * vc[1] +
    #     # 0 * vc[2] +
    #     7571/16695 * vc[3] +
    #     393/640 * vc[4] +
    #     -92097/339200 * vc[5] +
    #     187/2100 * vc[6] +
    #     1/40 * vc[7]
    # )
    #
    # @. vs = vc[1]

    vs
end
