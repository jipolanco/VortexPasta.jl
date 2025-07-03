using ..Constants: Infinity
using LinearAlgebra: ⋅

# Return the shortest separation r′ given the separation r = a - b between two points given a period L.
# In the periodic line, the shortest separation satisfies |r′| ≤ L / 2.
@inline function deperiodise_separation(r::Real, L::Real, Lhalf::Real)
    while r > Lhalf
        r -= L
    end
    while r < -Lhalf
        r += L
    end
    # @assert abs(r) ≤ Lhalf
    r
end

# Infinite period -> nothing to deperiodise.
@inline deperiodise_separation(r::Real, ::Infinity, ::Infinity) = r

# We convert SVector to tuple to make sure that no heap allocations are performed.
@inline function deperiodise_separation(r⃗::Vec3, args...)::Vec3
    r⃗ₚ = deperiodise_separation(Tuple(r⃗), map(Tuple, args)...)
    oftype(r⃗, r⃗ₚ)
end

@inline deperiodise_separation(rs::Tuple, args...)::Tuple = map(deperiodise_separation, rs, args...)

"""
    find_min_distance(
        fx::AbstractFilament, fy::AbstractFilament, i::Int, j::Int;
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
        maxiter = 4, rtol = 1e-2,
    )

Determine the minimum distance between filament segments.

This function estimates the minimum distance between filament segments `fx[i:i+1]` and `fy[j:j+1]`
via an iterative (Newton–Raphson) method.

The filaments `fx` and `fy` may be the same filament.

# Returns

Returns a `NamedTuple` with the following fields:

- `ζx`, `ζy`: optimal locations in ``[0, 1]`` within each segment;

- `x⃗`, `y⃗`: optimal locations within each segment;

- `p⃗`: periodic offset (each component is a multiple of the domain period along that direction);

- `d⃗`: minimum distance vector, `d⃗ = x⃗ - y⃗ + p⃗`.

# Optional keyword arguments

- `periods`: the period of the spatial domain. This should be given if one wants to take
  into account periodic images of the filaments;

- `maxiter = 4`: maximum number of iterations of the Newton method;

- `rtol = 1e-2`: relative tolerance for stopping the iterations. By default it's not very
  small since in general we don't need to be very precise.

"""
function find_min_distance(
        fx::AbstractFilament{T}, fy::AbstractFilament{T},
        i::Int, j::Int;
        periods::NTuple{3, Real} = (Infinity(), Infinity(), Infinity()),
        maxiter = 4,
        rtol::T = T(1e-2),
    ) where {T}
    periods_half = map(L -> L / T(2), periods)
    tx, ty = knots(fx), knots(fy)
    Δts = SVector(tx[i + 1] - tx[i], ty[j + 1] - ty[j])

    ζs = SVector(T(0.5), T(0.5))  # initial guess: middle of each segment

    # Use Newton's method to find minimum distance.
    # We minimise the squared distance d² between the two segments.
    for _ ∈ 1:maxiter
        ζx, ζy = ζs
        d⃗ = deperiodise_separation(fx(i, ζx) - fy(j, ζy), periods, periods_half)
        X′ = fx(i, ζx, Derivative(1))
        Y′ = fy(j, ζy, Derivative(1))
        X″ = fx(i, ζx, Derivative(2))
        Y″ = fy(j, ζy, Derivative(2))

        # Note: there's a factor 2 missing from the gradient and the Hessian,
        # but it's ok because the factor 2 would cancel out anyway when doing H \ ∇d².
        ∇d² = SVector(d⃗ ⋅ X′, -(d⃗ ⋅ Y′))

        # Hessian matrix
        Hxy = -(X′ ⋅ Y′)
        H = SMatrix{2, 2}(
            sum(abs2, X′) + (d⃗ ⋅ X″),
            Hxy, Hxy,
            sum(abs2, Y′) - (d⃗ ⋅ Y″)
        )

        # Note that these are negative increments (we want to do ζ_new = ζ - dζ).
        dts = H \ ∇d²
        dζs = dts ./ Δts  # normalise by knot increments (≈ segment lengths)

        # Make sure we stay within the segments (this is basically a constrained optimisation problem...).
        ζs_prev = ζs
        ζs = clamp.(ζs - dζs, 0, 1)
        dζs_actual = ζs - ζs_prev

        # Check if we converged to a minimum distance (we don't need to be very precise).
        maximum(abs, dζs_actual) < rtol && break
    end

    ζx, ζy = ζs
    x⃗ = fx(i, ζx)
    y⃗ = fy(j, ζy)
    r⃗ = x⃗ - y⃗
    d⃗ = deperiodise_separation(r⃗, periods, periods_half)
    p⃗ = d⃗ - r⃗  # periodic offset

    (; ζx, ζy, d⃗, x⃗, y⃗, p⃗,)
end
