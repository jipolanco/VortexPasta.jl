export NaiveShortRangeBackend

"""
    NaiveShortRangeBackend <: ShortRangeBackend

Naive computation of short-range interactions.

For now this is the only available short-range backend. In the future a more
efficient alternative backend will be added using cell lists.
"""
struct NaiveShortRangeBackend <: ShortRangeBackend end

struct NaiveShortRangeCache{
        Params <: ParamsShortRange,
        Timer <: TimerOutput,
    } <: ShortRangeCache
    params :: Params
    to     :: Timer
end

function init_cache_short(
        common::ParamsCommon{T}, params::ParamsShortRange{<:NaiveShortRangeBackend},
        to::TimerOutput,
    ) where {T}
    NaiveShortRangeCache(params, to)
end

# TODO split function into "elemental" parts that can be shared with another backend?
function short_range_velocity(
        g::F,  # this allows to set a custom kernel function g(α * r)
        cache::NaiveShortRangeCache,
        x⃗::Vec3,
        f::AbstractFilament,
        js::AbstractUnitRange = eachindex(segments(f)),
    ) where {F <: Function}
    (; params,) = cache
    (; common, quad, rcut,) = params
    (; Ls, Γ, α,) = common
    Lhs = map(L -> L / 2, Ls)  # half periods
    v⃗ = zero(x⃗)
    ts = knots(f)
    r²_cut = rcut * rcut
    Xs = nodes(f)
    t₊ = ts[first(js)]
    is_outside_range⁺ = let
        r⃗₊ = deperiodise_separation(x⃗ - Xs[first(js)], Ls, Lhs)
        r²₊ = sum(abs2, r⃗₊)
        r²₊ > r²_cut
    end
    ζs, ws = quadrature(quad)
    @inbounds for j ∈ js
        # Integrate over segment [j, j + 1]
        t₋, t₊ = t₊, ts[j + 1]
        r⃗₊ = deperiodise_separation(x⃗ - Xs[j + 1], Ls, Lhs)
        r²₊ = sum(abs2, r⃗₊)
        is_outside_range⁻, is_outside_range⁺ = is_outside_range⁺, r²₊ > r²_cut
        if is_outside_range⁻ && is_outside_range⁺
            # Skip this segment if its two extremities are too far from x⃗.
            continue
        end
        Δt = t₊ - t₋  # parameter increment
        for (ζ, w) ∈ zip(ζs, ws)
            X = f(j, ζ)
            Ẋ = f(j, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
            r⃗ = deperiodise_separation(x⃗ - X, Ls, Lhs)
            r² = sum(abs2, r⃗)
            r = sqrt(r²)
            v⃗ = v⃗ + (Δt * w * g(α * r) / r^3) * (Ẋ × r⃗)
        end
    end
    (Γ / 4π) * v⃗
end