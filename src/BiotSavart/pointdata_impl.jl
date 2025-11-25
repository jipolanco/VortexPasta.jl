"""
    set_num_points!(data::PointData, Np::Integer, quad::StaticSizeQuadrature)

Set the total number of non-uniform points that the cache must hold.

This will reallocate space to make all points fit in the cache. It will also reset the
contributions of previously-added charges.
"""
function set_num_points!(data::PointData, Np, quad::StaticSizeQuadrature)
    @assert Np < typemax(eltype(data.node_idx_prev))  # check that index type (e.g. Int32) is large enough
    Nq = Np * length(quad)  # number of quadrature nodes
    resize_no_copy!(data.nodes, Np)
    resize_no_copy!(data.node_idx_prev, Np)
    foreach(v -> resize_no_copy!(v, Np), data.derivatives_on_nodes)
    foreach(v -> resize_no_copy!(v, Np), data.subsegment_lengths)
    resize_no_copy!(data.points, Nq)
    resize_no_copy!(data.charges, Nq)
    data
end

# Total number of independent nodes among all filaments
_count_nodes(fs::AbstractVector{<:AbstractFilament}) = sum(f -> length(nodes(f)), fs; init = 0)

"""
    add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, params::ParamsBiotSavart)
    add_point_charges!(cache::LongRangeCache, fs::AbstractVector{<:AbstractFilament})

Add vector charges at multiple non-uniform locations.

This can be used for both short-range and long-range computations.

In periodic domains, each point is folded into the main periodic cell (in ``[0, L]^3``).
This can be accounted for to speed-up operations done later with point data, in particular in
short-range computations when determining the minimal distance between two points in the
periodic lattice.

In the case of long-range computations, this must be done before type-1 NUFFTs, to transform
from non-uniform data in physical space to uniform data in Fourier space. It must be called
before [`compute_vorticity_fourier!`](@ref).
"""
function add_point_charges!(data::PointData, fs::AbstractVector{<:AbstractFilament}, params::ParamsBiotSavart)
    Np = _count_nodes(fs)
    set_num_points!(data, Np, params.quad)
    chunks = FilamentChunkIterator(fs)  # defined in shortrange/shortrange.jl
    @sync for chunk in chunks
        isempty(chunk) && continue  # don't spawn a task if it will do no work
        Threads.@spawn let
            prev_indices = firstindex(fs):(first(chunk) - 1)  # filament indices given to all previous chunks
            n = _count_nodes(view(fs, prev_indices))  # we will start writing at index n + 1
            for i in chunk
                n = _add_point_charges!(data, fs[i], n, params, params.quad)
            end
        end
    end
    nothing
end

@inline function _compute_subsegment_lengths(f::AbstractFilament, i::Integer, quad, lims)
    δ⁻ = integrate(f, i - 1, quad; limits = lims[1]) do f, j, ζ
        @inline
        norm(f(j, ζ, Derivative(1)))
    end
    δ⁺ = integrate(f, i, quad; limits = lims[2]) do f, j, ζ
        @inline
        norm(f(j, ζ, Derivative(1)))
    end
    δ⁻, δ⁺
end

function _add_point_charges!(data::PointData, f::ClosedFilament, n::Int, params::ParamsBiotSavart, quad::StaticSizeQuadrature)
    (; Ls) = params
    (; lia_segment_fraction) = params.shortrange
    @assert eachindex(data.nodes) == eachindex(data.node_idx_prev) == eachindex(data.derivatives_on_nodes[1]) == eachindex(data.subsegment_lengths[1])
    @assert eachindex(data.points) == eachindex(data.charges)
    m = length(quad) * n  # current index in points/charges vectors
    nfirst = n + 1
    nlast = n + length(segments(f))
    mlast = length(quad) * nlast
    checkbounds(data.nodes, nlast)
    checkbounds(data.points, mlast)
    ζs, ws = quadrature(quad)
    ts = knots(f)
    subsegment_lims = lia_integration_limits(lia_segment_fraction)
    @inbounds for i in eachindex(nodes(f))
        n += 1
        data.nodes[n] = Filaments.fold_coordinates_periodic(f[i], Ls)
        data.node_idx_prev[n] = n - 1
        data.derivatives_on_nodes[1][n] = f[i, Derivative(1)]
        data.derivatives_on_nodes[2][n] = f[i, Derivative(2)]
        subsegs = _compute_subsegment_lengths(f, i, quad, subsegment_lims)
        data.subsegment_lengths[1][n] = subsegs[1]
        data.subsegment_lengths[2][n] = subsegs[2]
        Δt = ts[i + 1] - ts[i]
        for (ζ, w) ∈ zip(ζs, ws)
            s⃗ = f(i, ζ)
            s⃗′ = f(i, ζ, Derivative(1))  # = ∂f/∂t (w.r.t. filament parametrisation / knots)
            # Note: the vortex circulation Γ is included in the Ewald operator and
            # doesn't need to be included here.
            q = w * Δt
            _add_pointcharge!(data, s⃗, q * s⃗′, Ls, m += 1)
        end
    end
    @inbounds data.node_idx_prev[nfirst] = nlast  # account for periodic wrapping (closed filaments)
    @assert n == nlast
    @assert m == mlast
    n
end

function _add_point_charges!(data::PointData, f::ClosedFilament, n::Int, params::ParamsBiotSavart, ::NoQuadrature)
    (; Ls) = params
    (; lia_segment_fraction) = params.shortrange
    @assert eachindex(data.nodes) == eachindex(data.points) == eachindex(data.charges) ==
        eachindex(data.derivatives_on_nodes[1]) == eachindex(data.subsegment_lengths[1])
    nfirst = n + 1
    nlast = n + length(segments(f))
    checkbounds(data.nodes, nlast)
    checkbounds(data.points, nlast)
    γ = something(lia_segment_fraction, true)  # true in the sense of 1 (if segment_fraction == nothing)
    @inbounds for i in eachindex(nodes(f))
        n += 1
        data.nodes[n] = Filaments.fold_coordinates_periodic(f[i], Ls)
        data.node_idx_prev[n] = n - 1
        data.derivatives_on_nodes[1][n] = f[i, Derivative(1)]
        data.derivatives_on_nodes[2][n] = f[i, Derivative(2)]
        data.subsegment_lengths[1][n] = γ * norm(f[i] - f[i - 1])
        data.subsegment_lengths[2][n] = γ * norm(f[i + 1] - f[i])
        s⃗ = (f[i] + f[i + 1]) ./ 2   # segment midpoint as quadrature node
        s⃗′_dt = f[i + 1] - f[i]
        _add_pointcharge!(data, s⃗, s⃗′_dt, Ls, n)
    end
    @inbounds data.node_idx_prev[nfirst] = nlast  # account for periodic wrapping (closed filaments)
    @assert n == nlast
    n
end

@inline function _add_pointcharge!(data::PointData, X::Vec3, Q::Vec3, Ls::NTuple{3}, i::Int)
    @inbounds data.points[i] = Filaments.fold_coordinates_periodic(X, Ls)
    @inbounds data.charges[i] = Q
    data
end

